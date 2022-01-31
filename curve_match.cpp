#include "curve_match.h"
#include "score.h"

#include <QDateTime>
#include <QElapsedTimer>
#include <QTextStream>
#include <QStringList>
#include <QFileInfo>
#include <cmath>
#include <iostream>
#include <math.h>
#include <string.h>
#include <stdio.h> // for rename()

#include <sys/time.h>
#include <sys/resource.h>

#include "functions.h"

#define PARAMS_IMPL
#include "params.h"

#include "kb_distort.h"

#define BUILD_TS (char*) (__DATE__ " " __TIME__)

#define MISC_ACCT(word, coef_name, coef_value, value) { if (abs(value) > 1E-5) { DBG("     [*] MISC %s %s %f %f", (word), (coef_name), (float) (coef_value), (float) (value)) } }

using namespace std;

double getCPUTime() {
  struct rusage usage;
  if (getrusage(RUSAGE_SELF, &usage)) { return 0; } // ignored
  return (double) usage.ru_utime.tv_sec + ((double) usage.ru_utime.tv_usec) / 1000000.0;
}

void curve_smooth(QList<CurvePoint> &curve) {
  int l = curve.size();
  for(int i = 2 ; i < l - 2; i ++) {
    curve[i].smoothx = (-3 * curve[i - 2].x + 12 * curve[i - 1].x + 17 * curve[i].x + 12 * curve[i + 1].x - 3 * curve[i + 2].x)/35;
    curve[i].smoothy = (-3 * curve[i - 2].y + 12 * curve[i - 1].y + 17 * curve[i].y + 12 * curve[i + 1].y - 3 * curve[i + 2].y)/35;
  }
  for(int i = 0 ; i < 2; i ++) {
    curve[i].smoothx = curve[i].x;
    curve[i].smoothy = curve[i].y;
    curve[l - 1 - i].smoothx = curve[l - 1 - i].x;
    curve[l - 1 - i].smoothy = curve[l - 1 - i].y;
  }
}

/* --- main class for curve matching ---*/
CurveMatch::CurveMatch() : keyShift(&params) {
  loaded = false;
  params = default_params;
  debug = false;
  done = false;
  kb_preprocess = true;
  id = -1;
  screen_x = screen_y = 0;
  pixels_x = pixels_y = 0;
}

bool CurveMatch::curvePreprocess1(int curve_id) {
  /* curve preprocessing that can be evaluated incrementally :
     - evaluate turn rate
     - find sharp turns
     - normal vector calculation */

  bool on_hold = false;

  QList<CurvePoint> oneCurve = QList<CurvePoint>();
  QList<int> idxmap = QList<int>();

  bool end_flag = false;
  CurvePoint lastPoint(Point(-1, -1), curve_id, -1);
  for(int i = 0; i < curve.size(); i++) {
    if (curve[i].curve_id != curve_id) { continue; }
    if (curve[i].end_marker) {
      end_flag = true;
      break;
    }

    // point deduplication (workaround for bad data in old replay logs)
    if (curve[i].x == lastPoint.x && curve[i].y == lastPoint.y) { continue; }

    oneCurve.append(curve[i]);
    idxmap.append(i);

    lastPoint = curve[i];
  }
  int l = oneCurve.size();

  if (l < 2) { return false; }

  // apply smoothing to this curve ("smoothed" values will be in .smoothx & .smoothy attributes)
  if (params.smooth) {
    curve_smooth(oneCurve);
  } else {
    for(int i = 0; i < l; i ++) {
      oneCurve[i].smoothx = oneCurve[i].x;
      oneCurve[i].smoothy = oneCurve[i].y;
    }
  }

  for (int i = 1; i < l - 1; i ++) {
    oneCurve[i].turn_angle = (int) int(angle(oneCurve[i].smoothx - oneCurve[i-1].smoothx,
					     oneCurve[i].smoothy - oneCurve[i-1].smoothy,
					     oneCurve[i+1].smoothx - oneCurve[i].smoothx,
					     oneCurve[i+1].smoothy - oneCurve[i].smoothy) * 180 / M_PI + 0.5); // degrees


  }
  // avoid some side effects on curve_score (because there is often a delay before the second point)
  if (l > 1) {
    oneCurve[0].turn_angle = 0;
    oneCurve[l-1].turn_angle = 0;

    for (int i = 1; i < l - 1 ; i ++) {
      int t1 = oneCurve[i-1].turn_angle;
      int t2 = oneCurve[i].turn_angle;
      int t3 = oneCurve[i+1].turn_angle;

      if (abs(t2) > 160 && t1 * t2 < 0 && t2 * t3 < 0) {
	oneCurve[i].turn_angle += 360 * ((t2 < 0) - (t2 > 0));
      }
    }

    for (int i = 1; i < l - 1 ; i ++) {
      oneCurve[i].turn_smooth = int(0.5 * oneCurve[i].turn_angle + 0.25 * oneCurve[i-1].turn_angle + 0.25 * oneCurve[i+1].turn_angle);
    }
    oneCurve[0].turn_smooth = oneCurve[1].turn_angle / 4;
    oneCurve[l-1].turn_smooth = oneCurve[l-2].turn_angle / 4;
  }

  for(int i = 0 ; i < l; i ++) {
    oneCurve[i].sharp_turn = 0;
    oneCurve[i].flags = 0;
  }

  if (l >= 8) { // do not process small curves, probably a simple 2-letter words

    /* rotation / turning points */
    int sharp_turn_index = -1;
    int last_total_turn = -1;
    int last_turn_index = -100;
    int range = 1;
    int tips_margin = 2;
    for(int i = tips_margin ; i < l - tips_margin; i ++) {
      float total = 0;
      float t_index = 0;
      for(int j = i - range; j <= i + range; j ++) {
	total += oneCurve[j].turn_angle;
	t_index += oneCurve[j].turn_angle * j;
      }

      if ((abs(total) < last_total_turn || i == l - tips_margin - 1) && last_total_turn > params.turn_threshold) {
	if (sharp_turn_index >= 2 && sharp_turn_index < l - 2) {

	  for(int j = i - range; j <= i + range; j ++) {
	    if (abs(oneCurve[j].turn_angle) > params.turn_threshold2) {
	      sharp_turn_index = j;
	    }
	  }

	  int st_value = 1 + (last_total_turn > params.turn_threshold2 ||
			      abs(oneCurve[sharp_turn_index].turn_angle) > params.turn_threshold3);

	  if (st_value == 1 && abs(oneCurve[sharp_turn_index].turn_smooth) < params.turn_threshold_st6) {
	    st_value = 6;
	  }

	  // the code below does not handle ST=6 points properly
	  // cf. (unfinished) example in commit b1d9c03d72af34ff946c097f57abc8e54f254356
	  int diff = sharp_turn_index - last_turn_index;
	  if (diff <= 1) {
	    sharp_turn_index = -1;
	  } else if (diff == 2) {
	    sharp_turn_index --;
	    oneCurve[sharp_turn_index - 1].sharp_turn = 0;
	  } else if (diff <= params.turn_min_gap) {
	    bool remove_old = false, remove_new = false;
	    int old_value = oneCurve[last_turn_index].sharp_turn;

	    if (old_value < st_value) { remove_old = true; }
	    else if (old_value > st_value) { remove_new = true; }
	    else if (st_value == 1) {
	      remove_new = (abs(oneCurve[sharp_turn_index].turn_smooth) < abs(oneCurve[last_turn_index].turn_smooth));
	      remove_old = ! remove_new;
	    }

	    DBG("Sharp-turns too close: %d[%d] -> %d[%d] ---> remove_old=%d remove_new=%d",
		last_turn_index, old_value, sharp_turn_index, st_value, remove_old, remove_new);

	    if (remove_old) { oneCurve[last_turn_index].sharp_turn = 0; }
	    if (remove_new) { sharp_turn_index = -1; }
	  }

	  if (sharp_turn_index >= 0) {
	    oneCurve[sharp_turn_index].sharp_turn = st_value;

	    last_turn_index = sharp_turn_index;

	    DBG("Special point[%d]=%d (last_total_turn=%d, turn_angle=%d:%d)",
		sharp_turn_index, oneCurve[sharp_turn_index].sharp_turn,
		last_total_turn, oneCurve[sharp_turn_index].turn_angle, oneCurve[sharp_turn_index].turn_smooth);

	    sharp_turn_index = -1;
	  }
	}
      }

      if (abs(total) > params.turn_threshold) {
	sharp_turn_index = int(0.5 + abs(t_index / total));
      }

      last_total_turn = abs(total);
    }

    // for some missed obvious sharp turn (it sometime occurs near curve tips)
    for(int i = 2 ; i < l - 2; i ++) {
      if (oneCurve[i].sharp_turn == 0 && oneCurve[i - 1].sharp_turn == 0 && oneCurve[i + 1].sharp_turn == 0 &&
	  abs(oneCurve[i].turn_angle) > 120) {
	oneCurve[i].sharp_turn = 2;
      }
    }

  }

  /* --- acceleration computation --- */
  // timestamps smoothing
  int ts[l], len[l];
  len[0] = ts[0] = 0;
  for(int i = 1; i < l; i ++) {
    len[i] = len[i - 1] + distance(oneCurve[i - 1].smoothx, oneCurve[i - 1].smoothy, oneCurve[i].smoothx, oneCurve[i].smoothy) / scaling_ratio;
    ts[i] = oneCurve[i].t;
  }

  int savitsky_golay_coefs9[] = { -21, 14, 39, 54, 59, 54, 39, 14, -21 };

  for(int i = 5; i < l - 4; i ++) {
    if (! len[i - 4]) { continue; }

    int total = 0, total_coef = 0;
    for(int j = i - 4; j <= i + 4; j ++) {
      total += savitsky_golay_coefs9[j - i + 4] * oneCurve[j].t / len[j];
      total_coef += savitsky_golay_coefs9[j - i + 4];
    }
    ts[i] = len[i] * total / total_coef;
  }

  // speed evaluation
  float xspeed[l], yspeed[l];
  for(int i = 0; i < l - 1; i ++) {
    float dt = ts[i + 1] - ts[i];
    if (dt > 0) {
      xspeed[i] = (oneCurve[i + 1].smoothx - oneCurve[i].smoothx) / dt / scaling_ratio;
      yspeed[i] = (oneCurve[i + 1].smoothy - oneCurve[i].smoothy) / dt / scaling_ratio;
    } else if (i == 0) {
      xspeed[i] = yspeed[i] = 0;
    } else {
      xspeed[i] = xspeed[i - 1];
      yspeed[i] = yspeed[i - 1];
    }
  }
  xspeed[l - 1] = xspeed[l - 2];
  yspeed[l - 1] = yspeed[l - 2];

  // speedometer (used to find slow down points)
  int time_interval = params.speed_time_interval;
  for(int i = 0; i < l; i ++) {
    int i1 = i, i2 = i;
    if (! i1) { i1 = 1; } // avoid first point in case of device lag (which is very common)

    while (i1 > 1 && ts[i1] > ts[i] - time_interval) { i1 --; }
    while (i2 < l - 1 && ts[i2] < ts[i] + time_interval) { i2 ++; }
    float speed = 0;
    if (i2 > i1 && ts[i2] > ts[i1]) {
      speed = distance(oneCurve[i1].smoothx, oneCurve[i1].smoothy,
		       oneCurve[i2].smoothx, oneCurve[i2].smoothy) / (ts[i2] - ts[i1]) / scaling_ratio;
    } else if (i > 0) {
      speed = (float) oneCurve[i - 1].speed / 1000;
    }

    oneCurve[i].speed = 1000.0 * speed;
  }

  // speed smoothing
  float xsmooth[l], ysmooth[l];
  for(int i = 0; i < l - 1; i ++) {
    int sc = 1;
    float sx = 0, sy = 0;

    sx += 1 * xspeed[i]; sy += 1 * yspeed[i];
    bool fl = true, fr = true;
    for (int k = 1; k <= 4; k ++) { // @todo parameter
      if ((! fl) || (i - k < 0) || (oneCurve[i - k + 1].sharp_turn == 2)) {
	fl = false;
      } else {
	sx += xspeed[i - k]; sy += yspeed[i - k]; sc ++;
      }

      if ((! fr) || (i + k > l - 2) || (oneCurve[i + k].sharp_turn == 2)) {
	fr = false;
      } else {
	sx += xspeed[i + k]; sy += yspeed[i + k]; sc ++;
      }
    }

    xsmooth[i] = sx / sc;
    ysmooth[i] = sy / sc;
  }
  xsmooth[l - 1] = xsmooth[l - 2];
  ysmooth[l - 1] = ysmooth[l - 2];


  // acceleration evalution
  int accel[l];
  accel[0] = 0;
  for(int i = 1; i < l; i ++) {
    float dt = ts[i] - ts[i - 1];
    if (dt > 0) {
      float ax = (xsmooth[i] - xsmooth[i - 1]) / dt;
      float ay = (ysmooth[i] - ysmooth[i - 1]) / dt;
      float spd = distance(0, 0, xsmooth[i], ysmooth[i]);

      // normalize for speed twice (because acceleration is also "proportional" to speed for a given turn pattern)
      oneCurve[i].d2x = ax * 10000 / (spd * spd);
      oneCurve[i].d2y = ay * 10000 / (spd * spd);

      accel[i] = distance(0, 0, oneCurve[i].d2x, oneCurve[i].d2y);
      oneCurve[i].lac = (float) (oneCurve[i].d2x * ysmooth[i] - oneCurve[i].d2y * xsmooth[i]) / spd;
    } else {
      accel[i] = 0;
      oneCurve[i].lac = 0;
    }
  }

  // debug output (@todo remove)
  for(int i = 0; i <  l; i ++) {
    DBG("speed: %3d [%6d] [%6d] S=%5d/A=%5d %s (%6.2f,%6.2f) (%6.2f, %6.2f) ---> [%d,%d] = %d (%d)",
	i, oneCurve[i].t, ts[i], oneCurve[i].speed, accel[i], (oneCurve[i].dummy?"*":" "),
	xspeed[i], yspeed[i], xsmooth[i], ysmooth[i], oneCurve[i].d2x, oneCurve[i].d2y,
	accel[i], oneCurve[i].lac);
  }

  if (l >= 8) { // do not process small curves, probably a simple 2-letter words

    // alternate way of finding turning points
    int cur = 0;
    int total, start_index, max_turn, max_index;
    int threshold = params.atp_threshold;
    int min_angle1 = params.atp_min_angle1;
    int min_turn1 = params.atp_min_turn1;
    int max_pts = params.atp_max_pts;
    int tip_gap = 0;
    int opt_gap = params.atp_opt_gap;
    int excl_gap = params.atp_excl_gap;
    float pt61 = params.atp_pt_61;

    bool st_found = false;
    for(int i = tip_gap; i < l - tip_gap; i ++) {
      int turn = oneCurve[i].turn_smooth; // take care of jagged curves
      if (cur) {
	st_found |= (oneCurve[i].sharp_turn != 0);
	if (abs(turn) >= threshold && cur * turn > 0 && i < l - 1 - tip_gap) {
	  // turn continues
	  total += turn;
	  if (abs(turn) > abs(max_turn)) {
	    max_turn = turn;
	    max_index = i;
	  }
	} else {
	  // turn end
	  int end_index = i;
	  if (i == l - 1 - tip_gap) { total += turn; }

	  for(int j = max(max_index - excl_gap, tip_gap); j < min(max_index + excl_gap, l - tip_gap); j ++) {
	    if (oneCurve[j].sharp_turn) { st_found = true; }
	  }
	  if (max_index >= 2 && max_index < l + 2) {
	    DBG("New turn %d->%d total=%d max_index=%d max_turn=%d st_found=%d",
		start_index, end_index, total, max_index, abs(oneCurve[max_index].turn_smooth), st_found);
	    if (abs(total) > min_angle1 &&
		abs(oneCurve[max_index].turn_smooth) >= min_turn1 &&
		end_index - start_index <= max_pts &&
		! st_found) {

	      int value;
	      if (abs(oneCurve[max_index].turn_smooth) * pt61 > abs(oneCurve[max_index - 1].turn_smooth) &&
		  abs(oneCurve[max_index].turn_smooth) * pt61 > abs(oneCurve[max_index + 1].turn_smooth)) {
		value = 1;
	      } else {
		value = 6; // position of the matching point is not obvious
	      }
	      if (start_index <= opt_gap || end_index >= l - opt_gap) { value = 5; }
	      DBG("Special point[%d]=%d (try 2)", max_index, value);
	      oneCurve[max_index].sharp_turn = value;
	    }
	  }
	  cur = 0;
	}
      } else if (abs(turn) > threshold) {
	cur = (turn > 0) - (turn < 0);
	total = turn;
	start_index = i;
	max_index = i;
	max_turn = turn;
	st_found = (oneCurve[i].sharp_turn != 0);
      }
    }


    // add turning points based on acceleration
#define ACC(x) (accel[x]) /* abs(oneCurve[x].lac) */
    for(int i = 1; i < l - 1; i ++) {
      int threshold1 = params.accel_threshold1;
      int threshold2 = params.accel_threshold2;
      int gap = params.accel_gap;
      float ratio = params.accel_ratio;

      if (ACC(i) < threshold1) { continue; }
      if (oneCurve[i].sharp_turn) { continue; }

      if (ACC(i - 1) < ACC(i) / 3 || ACC(i + 1) < ACC(i) / 3) { continue; } // filter local "shaking"
      if (abs(oneCurve[i].lac) < ACC(i) / 3) { continue; } // filter acceleration with no radial component

      int min1 = 0, min2 = 0;
      bool ok = true;
      for (int j = 1; j <= gap; j ++) {
	if (i - j > 0) {
	  if (oneCurve[i - j].sharp_turn) { ok = false; break; }
	  if (ACC(i - j) > ACC(i)) { ok = false; break; }
	  if (ACC(i - j) < min1 || ! min1) { min1 = ACC(i - j); }
	}
	if (i + j < l) {
	  if (oneCurve[i + j].sharp_turn) { ok = false; break; }
	  if (ACC(i + j) > ACC(i)) { ok = false; break; }
	  if (ACC(i + j) < min2 || ! min2) { min2 = ACC(i + j); }
	}
      }

      if ((! ok) || (min1 > accel[i] * ratio) || (min2 > accel[i] * ratio)) { continue; }

      oneCurve[i].sharp_turn = (ACC(i) > threshold2)?1:5;
      DBG("Special point[%d]=%d (acceleration: %d - %d,%d)", i, oneCurve[i].sharp_turn, ACC(i), threshold1, threshold2);
    }

    // ugly workaround : adjust special points (1 & 3) that near missed an sharp angle turn
    for(int i = 2; i < l - 2; i ++) {
      if (oneCurve[i].sharp_turn != 1 && oneCurve[i].sharp_turn != 5) { continue; }
      for(int j = i - 1; j <= i + 1; j += 2) {
	if (abs(oneCurve[j].turn_angle > 50) and abs(oneCurve[j].turn_angle > 2 * oneCurve[i].turn_angle)) {
	  DBG("Special point moved: %d->%d (ST=%d)", i, j, oneCurve[i].sharp_turn);
	  oneCurve[j].sharp_turn = oneCurve[i].sharp_turn;
	  oneCurve[i].sharp_turn = 0;
	  break;
	}
      }
    }

    // catch some leftover obvious special points (partial slow down)
    int gap = params.accel_gap;
    for(int i = gap; i < l - gap; i ++) {
      if (oneCurve[i].sharp_turn) { continue; }
      bool ok = true;
      unsigned char m = 0;
      for(int j = 1; j <= gap; j ++) {
	if (oneCurve[i + j].speed < oneCurve[i].speed ||
	    oneCurve[i - j].speed < oneCurve[i].speed) { ok = false; break; }
	if (oneCurve[i + j].sharp_turn || oneCurve[i - j].sharp_turn) { ok = false; break; }

	if (oneCurve[i + j].speed > oneCurve[i].speed * 1.5) { m |= 1; }
	if (oneCurve[i - j].speed > oneCurve[i].speed * 1.5) { m |= 2; }
      }
      if (ok && m == 3 &&
	  abs(oneCurve[i].turn_smooth + oneCurve[i - 1].turn_smooth + oneCurve[i + 1].turn_smooth) > 30) {
	DBG("Partial slow down ST[%d] = 5", i);
	oneCurve[i].sharp_turn = 5;
      }
    }

    // add intermediate special points where the curve is the most far away from the theoretical segment
    // this is a bit redundant with above acceleration based part, but it stills find some new turns
    // (but it probably could be removed without causing any harm)
    int i0 = 0;
    for(int i = 0; i < l; i ++) {
      int st1 = oneCurve[i].sharp_turn;
      if (st1 == 1 || st1 == 2 || st1 == 6) {
	if (i0 && (i - i0) > params.max_turn_index_gap) {
	  int st0 = oneCurve[i0].sharp_turn;
	  if (st0 == 1 || st0 == 2 || st0 == 6) {
	    int index = 0;
	    float dmax = 0;
	    for(int j = i0 + 1; j < i; j ++) {
	      float d = dist_line_point(oneCurve[i0], oneCurve[i], oneCurve[j]);

	      // do not add a ST=5 point if there is there is already some
	      // special point between two ST=1 / ST=2 points
	      if (oneCurve[j].sharp_turn > 2) { dmax = 0; break; }

	      if (d > dmax) {
		dmax = d;
		index = j;
	      }
	    }
	    if (index > i0 + params.min_turn_index_gap &&
		index < i - params.min_turn_index_gap &&
		dmax > params.inter_pt_min_dist) {
	      oneCurve[index].sharp_turn = 5;
	      DBG("Special point[%d] = 5 (intermediate point [%d,%d], distance=%d)", index, i0, i, (int) dmax);
	    }
	  }
	}
	i0 = i;
      }
    }

    // slow down point search
    int maxd = params.max_turn_index_gap;
    for(int i = maxd / 2; i < l - maxd / 2; i ++) {
      if (oneCurve[i].sharp_turn) { continue; }

      int spd0 = oneCurve[i].speed;
      int ok = 0;
      int index[2] = { -1, -1 };
      for(int dir = -1; dir <= 1; dir += 2) {
	for(int j = 1; j <= maxd; j ++) {
	  if (i + dir * j < 0 || i + dir * j >= l) { continue; }

	  int spd = oneCurve[i + dir * j].speed;
	  if (spd < spd0 || oneCurve[i + dir * j].sharp_turn) { ok = 0; break; }

	  if ((spd > spd0 * params.slow_down_ratio) ||
	      (spd0 == 0 && spd > spd0)) { // case of 0 speed point (if user does move during speed calculation window)
	    ok |= (1 << (dir>0));
	    if (index[(dir + 1) >> 1] == -1) {
	      index[(dir + 1) >> 1] = i + dir * j;
	    }
	  }
	}
      }
      int delay = ts[index[1]] - ts[index[0]];
      if ((ok == 3) &&
	  (abs(oneCurve[i].turn_smooth) <= params.slow_down_max_turn)) {
	oneCurve[i].sharp_turn = 3;
	DBG("Special point[%d]=3 indexed=[%d,%d] delay=%d", i, index[0], index[1], delay);
      }
    }

    // special points adjustment based on speed (ST=1 only)
    /* @todo try again after smoothing :-)
    for(int i = 3; i < l - 3; i ++) {
      int st = oneCurve[i].sharp_turn;
      if (st == 1) {
	int speed = oneCurve[i].speed;
	int found = -1;
	for(int j = i - 3; j <= i + 3; j ++) {
	  if (oneCurve[j].speed < speed) {
	    found = j;
	    speed = oneCurve[j].speed;
	  }
	}
	if (found == i - 3 || found == i + 3) { continue; } // not a local minima

	// avoid moving a special point near another one
	if (found > 0) {
	  for (int j = max(found - 3, 0); j <= min(found + 3, l - 1); j ++) {
	    if (j != i && oneCurve[j].sharp_turn) { found = -1; break; }
	  }
	}

	// do not move from obvious turns
	if (found > 0 && abs(oneCurve[found].turn_angle) < abs(oneCurve[i].turn_angle) / 2) { found = -1; }

	if (found > 0) {
	  DBG("Special point adjustment [%d]: %d -> %d", st, i, found);
	  oneCurve[i].sharp_turn = 0;
	  oneCurve[found].sharp_turn = st;
	  if (found > i) { i = found + 1; }
	}
      }
    }
    */

    // compute "normal" vector for turns (really lame algorithm)
    for(int i = 2; i < l - 2; i ++) {
      if (oneCurve[i].sharp_turn) {
	int sharp_turn_index = i;

	int i1 = sharp_turn_index - 1;
	int i2 = sharp_turn_index + 1;
	float x1 = oneCurve[i1].smoothx - oneCurve[i1 - 1].smoothx;
	float y1 = oneCurve[i1].smoothy - oneCurve[i1 - 1].smoothy;
	float x2 = oneCurve[i2 + 1].smoothx - oneCurve[i2].smoothx;
	float y2 = oneCurve[i2 + 1].smoothy - oneCurve[i2].smoothy;
	float l1 = sqrt(x1 * x1 + y1 * y1);
	float l2 = sqrt(x2 * x2 + y2 * y2);
	oneCurve[sharp_turn_index].normalx = 100 * (x1 / l1 - x2 / l2); // integer vs. float snafu -> don't loose too much precision
	oneCurve[sharp_turn_index].normaly = 100 * (y1 / l1 - y2 / l2);
      }
    }

    /* detect loop hints */
    if (params.hints_master_switch) {
      int i0 = 0;
      int small_count = 0;

      for(int i = 0; i < l; i ++) {
	int turn = oneCurve[i].turn_smooth;
	if (i == l - 1 && ! end_flag && i0) {
	  DBG("[Hint-O] possible partial loop -> keeping curve on hold [%d:%d]", i0, i);
	  on_hold = true;
	}

	if (i < l - 1 && abs(turn) >= ((i0 && small_count < params.hint_o_max_small)?params.hint_o_turn_min_middle:params.hint_o_turn_min)) {
	  if ((! i0) || (i0 && oneCurve[i0].turn_smooth * turn < 0)) { i0 = i; small_count = 0; }

	  if (i0 && abs(turn) < params.hint_o_turn_min) {
	    small_count ++;
	  } else {
	    small_count = 0;
	  }

	} else if (i0) {
	  int i1 = i - 1;

	  // Rewind to last significant turn
	  while(i1 > i0 && abs(oneCurve[i1].turn_smooth) <= params.hint_o_turn_min) { i1 --; }

	  // is the loop matching the beginning of the curve ?
	  bool start_ok = true;
	  for(int j = 2; j < i0 ; j++) {
	    if (abs(oneCurve[j].turn_smooth) < params.hint_o_turn_min / 2 ||
		oneCurve[j].turn_smooth * oneCurve[i0].turn_smooth < 0) { start_ok = false; }
	  }
	  if (start_ok) { i0 = 0; }

	  // is the loop matching the tail of the curve ?
	  bool end_ok = true;
	  for(int j = i1 + 1; j < l - 2 ; j++) {
	    if (abs(oneCurve[j].turn_smooth) < params.hint_o_turn_min / 2 ||
		oneCurve[j].turn_smooth * oneCurve[i0].turn_smooth < 0) { end_ok = false; }
	  }
	  if (end_ok) { i1 = l - 1; }

	  int total = 0;
	  for(int i = i0; i <= i1; i ++) { total += oneCurve[i].turn_smooth; }

	  DBG("[Hint-O] Candidate loop [%d:%d] total=%d tips=%d:%d", i0, i1, total, start_ok, end_ok);

	  if (abs(total) > params.hint_o_total_min &&
	      (i1 - i0) >= params.hint_o_min_segments) {
	    int xavg = 0, yavg = 0;
	    for(int j = i0; j <= i1; j ++) {
	      xavg += oneCurve[j].x;
	      yavg += oneCurve[j].y;
	    }
	    xavg /= (i1 - i0 + 1);
	    yavg /= (i1 - i0 + 1);
	    int found = -1;
	    int min_dist = 0;
	    int max_dist = 0;
	    for(int j = i0; j <= i1; j ++) {
	      int dist = distance(xavg, yavg, oneCurve[j].x, oneCurve[j].y);
	      if (dist < min_dist || ! min_dist) {
		min_dist = dist;
		found = j;
	      }
	      if (dist > max_dist) { max_dist = dist; }
	    }

	    if (min_dist >= 0 && max_dist <= params.hint_o_max_radius / scaling_ratio) {

	      // remove ST=2 caused by the loop
	      for(int j = i0; j <= i1; j ++) {
		if (oneCurve[j].sharp_turn != 2 ||
		    oneCurve[j].turn_angle < params.hint_o_spare_st2_angle ||
		    abs(j - found) <= params.hint_o_spare_st2_gap) { // <- if this test is not discriminant enough, replace turns by ST=5
		  oneCurve[j].sharp_turn = 0;
		}
	      }

	      for(int j = i0; j <= i1; j ++) {
		oneCurve[j].flags |= FLAG_HINT_o;
	      }
	      oneCurve[found].sharp_turn = 2;
	      oneCurve[found].flags |= FLAG_HINT_O;

	      DBG("[Hint-O] Loop detected at curve index %d (error=%d, radius=%d, total_turn=%d, indexes=[%d:%d])",
		  found, min_dist, max_dist, total, i0, i1);

	    } else {
	      DBG("[Hint-O] bad distance: %d < %d < %d", min_dist, max_dist, params.hint_o_max_radius);
	    }
	  }
	  i0 = 0;
	}
      }
    } /* loop hint */

    /* detect potential "small hop" patterns (useful for horizontal stokes which are very hard to discriminate) */
    if (params.hints_master_switch) {
      int range = params.hint_v_range;
      int maxgap = params.hint_v_maxgap;
      for(int i = maxgap; i < l - maxgap; i ++) {
	if (! oneCurve[i].sharp_turn) { continue; }
	if (oneCurve[i].flags) { continue; }

	// @todo remove this test:
	if ((oneCurve[i + maxgap].x - oneCurve[i].x) * (oneCurve[i - maxgap].x - oneCurve[i].x) >= 0) { continue; }

	if (oneCurve[i + maxgap].y > oneCurve[i].y || oneCurve[i - maxgap].y > oneCurve[i].y) { continue; }

	int turn0 = oneCurve[i].turn_smooth;
	int turn[] = { 0, turn0, 0 };
	int failed = false;
	for(int dir = -1; dir <= 1; dir += 2) {
	  int k = i + dir;
	  bool changed = false;
	  while(k > 0 && k < l - 1 && abs(k - i) <= range) {
	    if (oneCurve[k].flags) { failed = true; }
	    int this_turn = oneCurve[k].turn_smooth;
	    if (this_turn * turn0 > 0 && ! changed) {
	      turn[1] += this_turn;
	    } else {
		turn[dir + 1] += this_turn;
		changed = true;
	    }
	    k += dir;
	  }
	}

	bool ok = (abs(turn[0]) > params.hint_v_minturn2 &&
		   abs(turn[1]) > params.hint_v_minturn &&
		   abs(turn[2]) > params.hint_v_minturn2 &&
		   ! failed);
	DBG("[Hint-V] Possible candidate at curve index %d (turn=[%d:%d:%d] other_hint=%d) -> %s",
	    i, turn[0], turn[1], turn[2], failed, ok?"*OK*":"(ignored)");
	if (ok) {
	  // note: we can not move the position of the special point because this is just a candidate
	  oneCurve[i].flags |= FLAG_HINT_V;
	}
      }
    } /* hint: small hop */

  } /* l >= 8 */


  if (end_flag) {
    oneCurve << EndMarker(curve_id);
  }

  // update aggregated curve (for logs, replay ...)
  for(int i = 0; i < idxmap.size(); i++) {
    curve[idxmap[i]] = oneCurve[i];
  }

  return on_hold;
}

void CurveMatch::curvePreprocess2() {
  /* curve preprocessing that can be deferred until final score calculation:
     - statistics
     - average speed
     - check if curve looks like a straight line (used later in scoring)
   */

  /* special points count */
  st.st_special = 0;
  for(int i = 0; i < curve.size(); i ++) {
    st.st_special += (curve[i].sharp_turn > 0);
  }

  /* check for straight line & compute speed */
  st.st_speed = 0;
  for (int i = 0; i < curve_count; i ++) {
    float sc1 = 0, sc2 = 0;
    int total = 0;

    QList<int> indexes;
    for (int j = 0; j < curve.size(); j ++) {
      if (curve[j].curve_id != i) { continue; }
      if (curve[j].end_marker) { continue; }
      indexes.append(j);
    }

    int i0 = indexes.first();
    int i1 = indexes.last();
    if (curve[i1].t > curve[i0].t) {
      st.st_speed += 1000.0 * curve[i1].length / (curve[i1].t - curve[i0].t);
    }

    int tip_soft = params.straight_tip;
    for(int k = tip_soft / 2; k < indexes.size() - tip_soft / 2; k ++) {
      int j = indexes[k];

      float coef = (k <= tip_soft || k >= indexes.size() - 1 - tip_soft)?0.5:1;
      int turn = coef * curve[j].turn_smooth;

      DBG("Straight score: Curve #%d:%d:%d %.3f -> %d", i, k, j, coef, turn);

      total += turn;
      sc1 = max(sc1, (float) abs(turn) / params.straight_max_turn);
    }
    sc2 = abs(total) * (0.35 + 0.65 * min(1, (float) indexes.size() / 250)) / params.straight_max_total; // @todo use parameters
    float straight = max(sc1, sc2);
    logdebug("Straight curve score: %.2f (%.2f, %.2f)", straight, sc1, sc2);
    quickCurves[i].straight = straight;
  }
}

void CurveMatch::setDebug(bool debug) {
  this -> debug = debug;
}
void CurveMatch::clearKeys() {
  keys.clear();
}

void CurveMatch::addKey(Key key) {
  scaling_ratio = 0; // we will need to re-compute this
  kb_preprocess = true;
  if (key.letter) {
    keys[key.label] = key;
    DBG("Key: '%s' %d,%d %d,%d (letter '%c')",
	QSTRING2PCHAR(key.label), key.x, key.y, key.width, key.height, key.letter);
  }
}

void CurveMatch::clearCurve() {
  curve.clear();
  done = false;
  memset(&st, 0, sizeof(st));
  curve_count = 0;
  memset(curve_started, 0, sizeof(curve_started));
}


/* Matching algorithm has been tuned with Jolla1 keyboard,
   So it will work with keyboards with different size (in pixels).
   Let's just find a scaling ratio to counteract this issue.
   ratios are > 1 if keyboard is bigger (pixel count) due to size and/or
   resolution. */

#define JOLLA1_KEY_DISTANCE 54
#define JOLLA1_DPI 245

#define MIN_DEVICE_WIDTH 40
#define MAX_DEVICE_WIDTH 80

#define MM_TO_INCHES 25.4

void CurveMatch::computeScalingRatio() {
  if (scaling_ratio != 0) { return; } // already computed or failed
  float total = 0;
  int count = 0;
  foreach(Key k1, keys) {
    float min_dist = 0;
    foreach(Key k2, keys) {
      if (k1.letter == k2.letter) { continue; }
      float d = distance(k1.x, k1.y, k2.x, k2.y);
      if (d < min_dist || ! min_dist) { min_dist = d; }
    }
    total += min_dist;
    count += 1;
  }

  if (! count) {
    // this may occur with some race condition
    scaling_ratio = -1;
    DBG("Scaling ratio compute error");
    return;
  }
  float avg_key_distance = total / count;

  float diagonal_inches = 0;
  int real_dpi = dpi;
  logdebug("real dpi: %.2f", real_dpi);
  if (real_dpi < 105) {
    // bogus screen information: probably obsolete SFOS1 device, old test of bad port
    // TODO remove this (require adaptation of all existing tests & detailed history)
    // try to guess based on screen size
    float diagonal_px = sqrt(pixels_x * pixels_x + pixels_y * pixels_y);
    logdebug("diagonal_px: %.2f", diagonal_px);
    if (diagonal_px > 1230 && diagonal_px < 2400) {
      // screen size may match Jolla C / Fairphone 2 / Nexus 5, let's assume it is a 5" screen
      // (this occurs at least on Nexus 5 with CM11 + SFOS 1.1.9)
      logdebug("Entered diagonal px");
      this->screen_size = diagonal_inches = 5;
      real_dpi = diagonal_px / diagonal_inches;
    } else {
      // fall back to Jolla1 and ignore all measurements (probably an old test) --> scaling_ratio = 1
      real_dpi = JOLLA1_DPI;
      avg_key_distance = JOLLA1_KEY_DISTANCE;
    }

  } else {
    // screen information available

    if (screen_x) {
      this->screen_size = diagonal_inches = sqrt(screen_x * screen_x + screen_y * screen_y) / MM_TO_INCHES;
      if (diagonal_inches >= 7.8 && diagonal_inches <= 7.9 && pixels_x == 1536) {
	// hardcoded mode for the tablet for people insisting on swiping on the tablet

	/* some dummy values:
           scaling_ratio = 1.33;
           dpi_ratio = size_ratio = sqrt(scaling_ratio);
	*/
      } else {
        /* Standard device: acceptance criteria based on screen width
           No more using diagonal because of different aspects ratios */
        if (screen_x < MIN_DEVICE_WIDTH || screen_x > MAX_DEVICE_WIDTH) {
          scaling_ratio = -1;
          dpi_ratio = size_ratio = 0;
        }
      }
    }
  }

  if (scaling_ratio == 0) {
    dpi_ratio = ((float) real_dpi) / JOLLA1_DPI;
    size_ratio = avg_key_distance / ((float) JOLLA1_KEY_DISTANCE) / dpi_ratio;

    /* Scale relatively to reference (Jolla1) size & resolution

       "params.scaling_size_pow" parameter value:
       0 : only compensate for resolution and assume all features scale linearly with keyboard size
       1 : compensate for resolution and size and assume all features size does not depend on keyboard physical size
       TODO: auto-tune this coefficient when tests case on other devices */

    scaling_ratio = dpi_ratio * pow(size_ratio, params.scaling_size_pow); // TODO temporary formula :-)

    scaling_ratio *= params.scaling_ratio_multiply;

    if (scaling_ratio > 0.95 && scaling_ratio < 1.05) {
      scaling_ratio = 1; // ignore a 5% change to avoid regression on existing tests
    }

  }

  params.glob_size_ratio = size_ratio;

  logdebug("PP: Average key distance: %.2f - DPI: %d - Diagonal: %.2f (%.1fx%.1f, %dx%d px)  %.2f (dpi: %.2f, size: %.2f)",
           avg_key_distance, real_dpi, diagonal_inches, screen_x, screen_y, pixels_x, pixels_y,
           scaling_ratio, dpi_ratio, size_ratio);

  if (params.scaling_ratio_override != 0) {
    scaling_ratio = params.scaling_ratio_override;
    logdebug("Scaling ratio override by configuration file: %.2f", scaling_ratio);
    return;
  }
}

void CurveMatch::setScreenInfo(int dpi, float screen_x, float screen_y) {
  this->dpi = dpi;
  this->screen_x = screen_x;
  this->screen_y = screen_y;
  this->scaling_ratio = 0;
  DBG("Screen DPI: %d - size: (%.2f, %.2f) (%dx%d px)", dpi, screen_x, screen_y, pixels_x, pixels_y);
}

void CurveMatch::setScreenSizePixels(int pixels_x, int pixels_y) {
  this->pixels_x = pixels_x;
  this->pixels_y = pixels_y;
}

void CurveMatch::addPoint(Point point, int curve_id, int timestamp) {
  QTime now = QTime::currentTime();

  if (curve_id < 0 || curve_id > MAX_CURVES) { return; }

  if (! scaling_ratio) {
    // compute scaling ratio
    computeScalingRatio();
  }

  if (kb_preprocess) {
    // we must apply keyboard biases before feeding the curve
    // in case of incremental processing
    kb_preprocess = false;
    if (params.thumb_correction) {
      float kb_scaling_ratio = dpi_ratio * pow(size_ratio, params.scaling_kb_size_pow); // TODO temporary formula :-)

      kb_distort(keys, params, kb_scaling_ratio);
    } else {
      kb_distort_cancel(keys);
    }
    keyShift.loadAndApply(keys);

    if (debug) {
      DBG("Keys adjustments:");
      QHashIterator<QString, Key> ki(keys);
      while (ki.hasNext()) {
	ki.next();
	Key key = ki.value();
	DBG("Key '%s' %d,%d -> %d,%d (letter='%c')",
	    QSTRING2PCHAR(key.label), key.x, key.y, key.corrected_x, key.corrected_y, key.letter);
      }
    }
  }

  if (curve.isEmpty()) {
    startTime = now;
    curve_length = 0;
  }

  int ts = (timestamp >= 0)?timestamp:startTime.msecsTo(now);

  if (curve.size() > 0) {
    if (ts < curve.last().t) {
      /* This is a workaround for a bug that produced bad test cases with
	 time going back. This should be eventually removed */
      for(int i = 0; i < curve.size(); i++) {
	curve[i].t = ts;
      }
    }

    /* Current implementation is dependant on good spacing between points
       (which is a design fault).
       Last SailfishOS has introduced some "micro-lag" (unnoticeable by
       user but it leads to "spaces" in the curve (... or maybe it's just
       some placebo effect).
       This is a simple workaround for this problem. */
    int last_index = -1;
    for(int j = curve.size() - 1; j >= 0; j --) {
      if (curve[j].curve_id == curve_id) { last_index = j; break; }
    }

    if (last_index >= 0) {
      CurvePoint lastPoint = curve[last_index];
      float dist = distancep(lastPoint, point);

      int max_length = params.max_segment_length;
      if (dist > max_length * scaling_ratio) {
	Point dp = point - lastPoint;
	int n = dist / (max_length * scaling_ratio);
	for(int i = 0; i < n; i ++) {
	  float coef = (i + 1.0) / (n + 1);
	  curve << CurvePoint(lastPoint + dp * coef, curve_id, lastPoint.t + coef * (ts - lastPoint.t), curve_length + coef * dist, true /* "dummy" point */);
	}
      }
      curve_length += dist;
    }

  }

  curve_started[curve_id] = true;
  curve_count = max(curve_count, curve_id + 1);

  curve << CurvePoint(point, curve_id, ts, curve_length);
}

void CurveMatch::endOneCurve(int curve_id) {
  if (curve_id < 0 || curve_id > MAX_CURVES) { return; }
  curve << EndMarker(curve_id);
  curve_started[curve_id] = false;
}



bool CurveMatch::loadTree(QString fileName) {
  /* load a .tre (word tree) file */

  keyShift.setDirectory(QFileInfo(fileName).path()); // by convention key-shift will use the same directory as .tre files

  if (loaded && fileName == this -> treeFile) { return true; }
  userDictionary.clear();

  bool status = wordtree.loadFromFile(fileName);
  loaded = status;
  this -> treeFile = fileName;
  if (fileName.isEmpty()) {
    logdebug("loadtree(-): %d", status);
    this -> userDictFile = QString();
    loaded = false;

  } else if (loaded) {
    QString uf = fileName;
    if (uf.endsWith(".tre")) { uf.remove(uf.length() - 4, 4); }
    uf += "-user.txt";
    this -> userDictFile = uf;
    loadUserDict();
    logdebug("loadTree(%s): %d", QSTRING2PCHAR(fileName), status);
  }
  return status;
}

void CurveMatch::log(QString txt) {
  if (! logFile.isNull()) {
    QFile file(logFile);
    if (file.open(QIODevice::Append)) {
      QTextStream out(&file);
      out << txt << "\n";
    }
    file.close();
  }
}

void CurveMatch::setLogFile(QString fileName) {
  logFile = fileName;

  log_setfile(fileName); // for debug logging
}

QList<ScenarioType> CurveMatch::getCandidates() {
  return candidates;
}

QList<ScenarioDto> CurveMatch::getCandidatesDto() {
  QList<ScenarioDto> result;
  foreach(ScenarioType s, candidates) {
    ScenarioDto dto(s.getNameRealLetters(), s.getWordListAsList().join(QString(",")),
		    s.getScore(), s.getClass(), s.getStar());
    result.append(dto);
  }
  return result;
}

void CurveMatch::scenarioFilter(QList<ScenarioType> &scenarios, float score_ratio, int min_size, int max_size, bool finished) {
  /* "skim" any scenario list based on number and/or score

     This method is awfully inefficient, but as the target is to use only the incremental
     implementation (class IncrementalMatch), it'll stay as is */

  float max_score = 0;

  qSort(scenarios.begin(), scenarios.end());  // <- 25% of the whole CPU usage, ouch!
  // (in incremental mode it's only used to sort candidates, so it's not a problem)

  foreach(ScenarioType s, scenarios) {
    float sc = s.getScore();
    if (sc > max_score) { max_score = sc; }
  }

  QHash<QString, int> dejavu;

  int i = 0;
  while (i < scenarios.size()) {
    float sc = scenarios[i].getScore();

    if (sc < max_score * score_ratio && scenarios.size() > min_size) {
      // remove scenarios with lowest scores
      st.st_skim ++;
      DBG("filtering(score): \"%s\" %s (%.3f/%.3f)", QSTRING2PCHAR(scenarios[i].getName()), QSTRING2PCHAR(scenarios[i].getId()), sc, max_score);
      scenarios.takeAt(i);

    } else if (finished || ! scenarios[i].forkLast()) {
      // remove scenarios with duplicate words (but not just after a scenario fork)

      QString name = scenarios[i].getName();
      if (dejavu.contains(name)) {
	int i0 = dejavu[name];
	scenarios.takeAt(i0);
	DBG("filtering(fork): \"%s\" %s (%.3f/%.3f)", QSTRING2PCHAR(scenarios[i0].getName()), QSTRING2PCHAR(scenarios[i0].getId()), sc, max_score);
	foreach(QString n, dejavu.keys()) {
	  if (dejavu[n] > i0) {
	    dejavu[n] --;
	  }
	}
	i --;
	dejavu.remove(name); // will be added again on next iteration
      } else {
	dejavu.insert(name, i);
	i++;
      }

    } else {
      i++;

    }

  }

  // enforce max scenarios count
  while (max_size > 1 && scenarios.size() > max_size) {
    st.st_skim ++;
    ScenarioType s = scenarios.takeAt(0);
    DBG("filtering(size): \"%s\" %s (%.3f/%.3f)", QSTRING2PCHAR(s.getName()), QSTRING2PCHAR(s.getId()), s.getScore(), max_score);
  }

}

bool CurveMatch::setCurves() {
  bool result = false;
  for (int i = 0; i < curve_count; i ++) {
    bool on_hold = curvePreprocess1(i);
    quickCurves[i].setCurve(curve, scaling_ratio, i, params.multi_dot_threshold);
    quickCurves[i].on_hold = on_hold;

    result |= on_hold;
  }
  quickCurves[curve_count].clearCurve();
  return result;
}

bool CurveMatch::match() {
  /* run the full "one-shot algorithm */
  scenarios.clear();
  candidates.clear();

  if (! loaded || ! keys.size() || ! curve.size()) { return false; }
  if (curve.size() < 3) { return false; }

  memset(& st, 0, sizeof(st));

  // change order for equal items: qSort(curve.begin(), curve.end()); // in multi mode we may lose point ordering

  quickKeys.setParams(&params);
  quickKeys.setKeys(keys, scaling_ratio);

  setCurves();
  curvePreprocess2();

  ScenarioType root = ScenarioType(&wordtree, &quickKeys, quickCurves, &params);
  root.setDebug(debug);
  scenarios.append(root);

  QElapsedTimer timer;
  timer.start();
  double start_cpu_time = getCPUTime();

  int count = 0;

  int n = 0;
  while (scenarios.size()) {
    QList<ScenarioType> new_scenarios;
    foreach(ScenarioType scenario, scenarios) {

      QList<ScenarioType> childs;
      scenario.nextKey(childs, st);
      foreach(ScenarioType child, childs) {
	if (child.isFinished()) {
	  if (child.postProcess(st)) {
	    DBG("New candidate: %s (score=%.3f)", QSTRING2PCHAR(child.getId()), child.getScore());
	    candidates.append(child);
	  } else {
	    DBG("Failed candidate: %s", QSTRING2PCHAR(child.getId()));
	  }
	} else {
	  new_scenarios.append(child);
	}
	count += 1;
      }

    }
    n += 1;
    scenarios = new_scenarios; // remember QT collections have intelligent copy-on-write (hope it works)

    if (debug) {
      foreach(ScenarioType scenario, scenarios) {
	DBG("LST[%d] %s [%.2f]", n, QSTRING2PCHAR(scenario.getId()), scenario.getScore());
      }
    }

    if (n >= 3) {
      // scenario filtering works differently in non-incremental mode
      // as we do not care about performance in this case, just arbitrary increase max_active_scenarios value
      scenarioFilter(scenarios, 0, 15, params.max_active_scenarios * 2, false); // @todo add to parameter list
      DBG("Depth: %d - Scenarios: %d - Candidates: %d", n, scenarios.size(), candidates.size());
    }
  }

  st.st_cputime = (int) (1000 * (getCPUTime() - start_cpu_time));
  st.st_time = (int) timer.elapsed();
  st.st_count = count;

  scenarioFilter(candidates, 0.5, 10, 3 * params.max_candidates, true); // @todo add to parameter list
  sortCandidates();
  scenarioFilter(candidates, 0.7, 10, params.max_candidates, true); // @todo add to parameter list

  storeKeyPos();


  logdebug("Candidates: %d (time=%d, nodes=%d, forks=%d, skim=%d, speed=%d, special=%d, cputime=%d, treefile=%s)",
	   candidates.size(), st.st_time, st.st_count, st.st_fork, st.st_skim, st.st_speed,
	   st.st_special, st.st_cputime, QSTRING2PCHAR(treeFile));

  done = true;

  return candidates.size() > 0;
}

void CurveMatch::sortCandidates() {
  QList <ScenarioType *> pcandidates = QList<ScenarioType *>();
  QListIterator<ScenarioType> it(candidates);
  while(it.hasNext()) {
    pcandidates.append((ScenarioType*) &(it.next()));
  }

  ScenarioType::sortCandidates(pcandidates, params, debug);
}

void CurveMatch::endCurve(int correlation_id) {
  for (int i = 0; i < MAX_CURVES; i ++) {
    if (curve_started[i]) {
      curve << EndMarker(i);
      DBG("Implicit curve end #%d", i);
    }
  }
  this -> id = correlation_id;
  log(QString("IN: ") + toString());
  if (! done) { match(); }
  log(QString("OUT: ") + resultToString());
}

void CurveMatch::setParameters(QString jsonStr) {
  QJsonDocument doc = QJsonDocument::fromJson(jsonStr.toUtf8());
  params = Params::fromJson(doc.object());
  kb_preprocess = true;
}

void CurveMatch::useDefaultParameters() {
  params = default_params;
}

QList<CurvePoint> CurveMatch::getCurve() {
  return curve;
}

void CurveMatch::toJson(QJsonObject &json) {
  json["id"] = id;

  // parameters
  QJsonObject json_params;
  params.toJson(json_params);
  json["params"] = json_params;

  // keys
  QJsonArray json_keys;
  QHash<QString, Key>::iterator i;
  for (i = keys.begin(); i != keys.end(); ++i) {
    Key k = i.value();
    QJsonObject json_key;
    k.toJson(json_key);
    json_keys.append(json_key);
  }
  json["keys"] = json_keys;
  json["key_hash"] = keyShift.getHash();

  // curve
  QJsonArray json_curve;
  foreach(CurvePoint p, curve) {
    QJsonObject json_point;
    p.toJson(json_point);
    json_curve.append(json_point);
  }
  json["curve"] = json_curve;

  // other
  json["treefile"] = treeFile;
  json["datetime"] = QDateTime::currentDateTime().toString(Qt::ISODate);

  // screen features
  computeScalingRatio(); // make sure it is up-to-date
  QJsonObject json_scaling;
  json_scaling["scaling_ratio"] = scaling_ratio;
  json_scaling["dpi_ratio"] = dpi_ratio;
  json_scaling["size_ratio"] = size_ratio;
  json_scaling["dpi"] = dpi;
  json_scaling["screen_size"] = screen_size;
  json_scaling["screen_x"] = screen_x;
  json_scaling["screen_y"] = screen_y;
  json_scaling["pixels_x"] = pixels_x;
  json_scaling["pixels_y"] = pixels_y;
  json["scaling"] = json_scaling;
}

void CurveMatch::fromJson(const QJsonObject &json) {
  QJsonValue input = json["input"];
  if (input != QJsonValue::Undefined) {
    fromJson(input.toObject());
    return;
  }

  // parameters
  params = Params::fromJson(json["params"].toObject());

  // keys
  QJsonArray json_keys = json["keys"].toArray();
  keys.clear();
  foreach(QJsonValue json_key, json_keys) {
    Key k = Key::fromJson(json_key.toObject());
    keys[k.label] = k;
  }
  scaling_ratio = 0;
  computeScalingRatio();

  // curve
  QJsonArray json_curve = json["curve"].toArray();
  curve.clear();
  foreach(QJsonValue json_point, json_curve) {
    QJsonObject json_object = json_point.toObject();
    CurvePoint p = CurvePoint::fromJson(json_object);
    if (! p.dummy) {
      curve.append(p);
    }
  }

  // scaling
  QJsonObject json_scaling = json["scaling"].toObject();
  dpi = json_scaling["dpi"].toDouble();
  screen_x = json_scaling["screen_x"].toDouble();
  screen_y = json_scaling["screen_y"].toDouble();
  pixels_x = json_scaling["pixels_x"].toDouble();
  pixels_y = json_scaling["pixels_y"].toDouble();
  scaling_ratio = 0; /* we will recompute it */
}

void CurveMatch::fromString(const QString &jsonStr) {
  QJsonDocument doc = QJsonDocument::fromJson(jsonStr.toUtf8());
  fromJson(doc.object());
}

QString CurveMatch::toString(bool indent) {
  QJsonObject json;
  toJson(json);
  QJsonDocument doc(json);
  return QString(doc.toJson(indent?QJsonDocument::Indented:QJsonDocument::Compact));
}

void CurveMatch::resultToJson(QJsonObject &json) {
  json["id"] = id;
  json["build"] = QString(BUILD_TS);

  QJsonObject json_input;
  toJson(json_input);
  json["input"] = json_input;

  QJsonArray json_candidates;
  QList<ScenarioType> candidates = getCandidates();
  qSort(candidates.begin(), candidates.end());
  foreach (ScenarioType c, candidates) {
    QJsonObject json_scenario;
    c.toJson(json_scenario);
    json_candidates.append(json_scenario);
  }
  json["candidates"] = json_candidates;
  json["ts"] = QDateTime::currentDateTime().toString(Qt::ISODate);

  QJsonObject json_stats;
  json_stats["time"] = st.st_time;
  json_stats["cputime"] = st.st_cputime;
  json_stats["count"] = st.st_count;
  json_stats["fork"] = st.st_fork;
  json_stats["skim"] = st.st_skim;
  json_stats["special"] = st.st_special;
  json_stats["speed"] = st.st_speed;
  json_stats["cache_hit"] = st.st_cache_hit;
  json_stats["cache_miss"] = st.st_cache_miss;
  json["stats"] = json_stats;

  QJsonObject json_params;
  default_params.toJson(json_params);
  json["default_params"] = json_params;

  QJsonArray json_st;
  for(int i = 0; i < curve_count; i ++) {
    json_st.append(quickCurves[i].straight);
  }
  json["straight"] = json_st;

  json["multi_count"] = curve_count;
}

QString CurveMatch::resultToString(bool indent) {
  QJsonObject json;
  resultToJson(json);
  QJsonDocument doc(json);
  return QString(doc.toJson(indent?QJsonDocument::Indented:QJsonDocument::Compact));
}

void CurveMatch::learn(QString word, int addValue, bool init) {
  QString letters = word2letter(word);

  DBG("CM-learn [%s]: %s += %d (init:%d)", QSTRING2PCHAR(letters), QSTRING2PCHAR(word), addValue, (int) init);

  if (letters.length() < 2 || word.length() < 2) { return; }
  if (! params.user_dict_learn) { return; }
  if (! loaded) { return; }

  QString payload_word = word;
  if (word == letters) {
    payload_word = QString("=");
  }

  // get user dictionary
  UserDictEntry entry;
  if (userDictionary.contains(word)) {
    entry = userDictionary[word];
  }

  // get data from in-memory tree
  QPair<void*, int> pl = wordtree.getPayload(QSTRING2PUCHAR(letters));

  // compute new node value for in-memory tree
  QString payload;
  bool found = false;
  if (pl.first) { // existing node
    payload = QString((const char*) pl.first);
    QStringList lst = payload.split(",");
    foreach(QString w, lst) {
      if (w == payload_word) {
	found = true; // we already know this word
      }
    }
    if (found) {
      payload = QString(); // don't update the tree
    } else {
      lst.append(payload_word);
      payload = lst.join(",");
    }
  } else { // new node
    payload = payload_word;
  }

  // do not learn new words if we already know them
  if (init && found) {
    // during initialization we try to add a word which is already in the tree
    // (may be caused by issues with storage file or updated dictionary)
    userDictionary.remove(word);
    return;
  }
  if (! init && ! userDictionary.contains(word) && found) {
    // during keyboard usage do not add known word unless they are already in user directory
    return;
  }

  // update in-memory tree
  int now = (int) time(0);
  if (addValue >= 0 && ! payload.isEmpty()) {
    if (! init) { logdebug("Learned new word: %s [%s]", QSTRING2PCHAR(word), QSTRING2PCHAR(letters)); }
    unsigned char *ptr = QSTRING2PUCHAR(payload);
    int len = strlen((char*) ptr) + 1;
    unsigned char pl_char[len];
    memmove(pl_char, ptr, len);
    wordtree.setPayload(QSTRING2PUCHAR(letters), pl_char, len);
    DBG("CM-learn: update tree %s -> '%s' payload=[%s]", QSTRING2PCHAR(letters), QSTRING2PCHAR(word), pl_char);
  }

  // update user directory
  float new_count;
  if (init) {
    new_count = entry.count;
  } else {
    new_count = entry.getUpdatedCount(now, params.user_dict_halflife) + addValue;
    if (new_count < 0) { new_count = 0; }

    userdict_dirty = true;
  }
  userDictionary[word] = UserDictEntry(letters, now, new_count);

  if (! init || debug) {
    logdebug("Learn: %s (%s) (init: %d, add: %d) %.4f->%.4f",
	     QSTRING2PCHAR(word), QSTRING2PCHAR(letters),
	     init, addValue,
	     entry.count, new_count);
  }
}

void CurveMatch::loadUserDict() {
  if (userDictFile.isEmpty()) { return; }

  userdict_dirty = false;
  userDictionary.clear();

  QFile file(userDictFile);
  if (! file.open(QIODevice::ReadOnly)) { return; }

  QTextStream in(&file);

  QString line;
  do {
    line = in.readLine();
    QTextStream stream(&line);
    QString word, letters;
    float count;
    int ts;

    stream >> word >> letters >> count >> ts;

    if (count && ts) {
      userDictionary[word] = UserDictEntry(letters, ts, count);
      learn(word, 0, true); // add the word to in memory tree dictionary
    }
  } while (!line.isNull());

  file.close();
  purgeUserDict();
}

void CurveMatch::saveUserDict() {
  if (! userdict_dirty) { return; }
  if (userDictFile.isEmpty()) { return; }
  if (userDictionary.isEmpty()) { return; }

  purgeUserDict();

  QFile file(userDictFile + ".tmp");
  if (! file.open(QIODevice::WriteOnly)) { return; }

  QTextStream out(&file);

  int now = (int) time(0);

  QHashIterator<QString, UserDictEntry> i(userDictionary);
  while (i.hasNext()) {
    i.next();
    UserDictEntry entry = i.value();
    QString word = i.key();

    float count = entry.getUpdatedCount(now, params.user_dict_halflife);

    if (count > params.user_dict_min_count && ! word.isEmpty() && ! entry.letters.isEmpty()) {
      out << word << " " << entry.letters << " " << entry.count << " " << entry.ts << endl;
    }
  }

  file.close();
  if (file.error()) { return; /* save failed */ }

  // QT rename can't do an atomic file replacement
  rename(QSTRING2PCHAR(userDictFile + ".tmp"), QSTRING2PCHAR(userDictFile));

  userdict_dirty = false;
}

float UserDictEntry::getUpdatedCount(int now, int days) const {
  if (! ts || ! count) { return 0; }
  if (now < ts) { return count; }

  return ((float) count) * exp(- (float) (now - ts) * log(2) / days / 86400);
}

static int userDictLessThan(QPair<QString, float> &e1, QPair<QString, float> &e2) {
  return e1.second < e2.second;
}


void CurveMatch::purgeUserDict() {
  if (userDictionary.size() <= params.user_dict_size) { return; }

  int now = time(0);

  userdict_dirty = true;

  QList<QPair<QString, float> > lst;
  QHashIterator<QString, UserDictEntry> i(userDictionary);
  while (i.hasNext()) {
    i.next();
    float score = i.value().getUpdatedCount(now, params.user_dict_halflife);
    lst.append(QPair<QString, float>(i.key(), score));
  }

  qSort(lst.begin(), lst.end(), userDictLessThan);

  for(int i = 0; i < lst.size() - params.user_dict_size; i ++) {
    userDictionary.remove(lst[i].first);
    logdebug("Learn/purge: %s", QSTRING2PCHAR(lst[i].first));
  }
}

void CurveMatch::dumpDict() {
  wordtree.dump();
}

char* CurveMatch::getPayload(unsigned char *letters) {
  return (char*) wordtree.getPayload(letters).first;
}

void CurveMatch::loadKeyPos() {
  // this is not done automatically when context is loaded from a JSON file
  // (this is used from CLI tool)
  keyShift.loadAndApply(keys);
}

void CurveMatch::storeKeyPos() {
  /* @TODO */
}

void CurveMatch::updateKeyPosForTest(QString /* expected */) {
  if (! params.key_shift_enable) { return; }
  /* broken for now
  DBG(KS_TAG "*test* Update for word %s", QSTRING2PCHAR(expected));

  int found = -1;
  for(int i = 0; i < candidates.size(); i ++) {
    QStringList words = candidates[i].getWordList().split(',');
    if (candidates[i].getName() == expected ||
	words.contains(expected)) {
      found = i;
    }
  }
  if (found == -1) {
    DBG(KS_TAG "*test* word not found !");
    return;
  }

  keyShift.lock();
  QList<QPair<unsigned char, Point> > match_points = candidates[found].get_key_error();

  for(int i = 0; i < match_points.size(); i ++) {
    // @TODO find the right key name as QString
    keyShift.update(match_points[i].first,
		    match_points[i].second.x,
		    match_points[i].second.y);
  }

  keyShift.save();
  DBG(KS_TAG "*test* update done !");
  */
}

void CurveMatch::saveKeyPos() {
  keyShift.save();
}

Params* CurveMatch::getParamsPtr() {
  return &params;
}
