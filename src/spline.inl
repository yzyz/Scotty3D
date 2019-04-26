// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
    const T& position0, const T& position1, const T& tangent0,
    const T& tangent1, double normalizedTime, int derivative) {
  // TODO (Animation) Task 1a
  double t = normalizedTime;
  double t2 = t * t;
  double t3 = t * t2;

  double h00 = 2 * t3 - 3 * t2 + 1;
  double h10 = t3 - 2 * t2 + t;
  double h01 = -2 * t3 + 3 * t2;
  double h11 = t3 - t2;
  return h00 * position0 + h10 * tangent0 + h01 * position1 + h11 * tangent1;
}

// Returns a state interpolated between the values directly before and after the
// given time.
template <class T>
inline T Spline<T>::evaluate(double time, int derivative) {
  // TODO (Animation) Task 1b
  if (knots.size() < 1)
    return T();
  else if (knots.size() == 1)
    return knots.begin()->second;
  else if (time <= knots.begin()->first)
    return knots.begin()->second;
  else if (time >= knots.rbegin()->first)
    return knots.rbegin()->second;
  else {
    double t0, t1, t2, t3;
    T p0, p1, p2, p3;

    auto it2 = knots.upper_bound(time);
    t2 = it2->first;
    p2 = it2->second;

    auto it1 = prev(it2);
    t1 = it1->first;
    p1 = it1->second;

    if (next(it2) == knots.end()) {
      t3 = t2 + t2 - t1;
      p3 = p2 + p2 - p1;
    }
    else {
      auto it3 = next(it2);
      t3 = it3->first;
      p3 = it3->second;
    }

    if (it1 == knots.begin()) {
      t0 = t1 - t2 + t1;
      p0 = p1 - p2 + p1;
    }
    else {
      auto it0 = prev(it1);
      t0 = it0->first;
      p0 = it0->second;
    }

    T m1 = (p2 - p0) / (t2 - t0) * (t2 - t1);
    T m2 = (p3 - p1) / (t3 - t1) * (t2 - t1);
    double t = (time - t1) / (t2 - t1);

    return cubicSplineUnitInterval(p1, p2, m1, m2, t);
  }
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance) {
  // Empty maps have no knots.
  if (knots.size() < 1) {
    return false;
  }

  // Look up the first element > or = to time.
  typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
  typename std::map<double, T>::iterator t1_iter;
  t1_iter = t2_iter;
  t1_iter--;

  if (t2_iter == knots.end()) {
    t2_iter = t1_iter;
  }

  // Handle tolerance bounds,
  // because we are working with floating point numbers.
  double t1 = (*t1_iter).first;
  double t2 = (*t2_iter).first;

  double d1 = fabs(t1 - time);
  double d2 = fabs(t2 - time);

  if (d1 < tolerance && d1 < d2) {
    knots.erase(t1_iter);
    return true;
  }

  if (d2 < tolerance && d2 < d1) {
    knots.erase(t2_iter);
    return t2;
  }

  return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue(double time, T value) {
  knots[time] = value;
}

template <class T>
inline T Spline<T>::operator()(double time) {
  return evaluate(time);
}
