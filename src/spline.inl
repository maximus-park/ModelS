// Given a time between 0 and 1, evaluates a cubic polynomial with
// the given endpoint and tangent values at the beginning (0) and
// end (1) of the interval.  Optionally, one can request a derivative
// of the spline (0=no derivative, 1=first derivative, 2=2nd derivative).
template <class T>
inline T Spline<T>::cubicSplineUnitInterval(
      const T& position0,
      const T& position1,
      const T& tangent0,
      const T& tangent1,
      double normalizedTime,
      int derivative )
{
   double t1 = normalizedTime;
   double t2 = normalizedTime * t1;
   double t3 = normalizedTime * t2;
   double h00, h10, h01, h11;
   
   switch (derivative) {
     default: // falling through
     case 0:
       h00 = 2.*t3 - 3.*t2 + 1.; h10 = t3 - 2.*t2 + t1;
       h01 = -2.*t3 + 3.*t2;     h11 = t3 - t2;
       break;
     case 1:
       h00 = 6.*t2 - 6.*t1;      h10 = 3.*t2 - 4.*t1 + 1.;
       h01 = -6.*t2 + 6.*t1;     h11 = 3.*t2 - 2.*t1; 
       break;
     case 2:
       h00 = 12.*t1 - 6.;        h10 = 6.*t1 - 4.;
       h01 = -12.*t1 + 6.;       h11 = 6.*t1 - 2.; 
       break;
   }
   return position0*h00 + tangent0*h10 + position1*h01 + tangent1*h11;
}
            
// Returns a state interpolated between the values directly before and after the given time.
template <class T>
inline T Spline<T>::evaluate( double time, int derivative )
{
   double t0, t1, t2, t3;
   T p0, p1, p2, p3;
   T m1, m2;
   KnotCIter prev, next, p_prev, n_next;
  
   // If there are no knots at all in the spline, interpolation
   // should return the default value for the interpolated type.
   if (knots.size() < 1) return T();
   
   // If there is only one knot in the spline,
   // interpolation should always return the value of
   // that knot (independent of the time). Derivatives should
   // be zero.
   if (knots.size() == 1) {
       if (derivative != 0) return T();
       else return knots.begin()->second;
   }
   // Get iterator to the first and last knot
   KnotCIter knot_first = knots.begin(),
	     knot_last  = std::prev(knots.end());
   
   // Clamp query time inside the spline domain
   // Find p1 p2 iterator
   if (time < knot_first->first) {
       time = knot_first->first;
       prev = knot_first;
       next = std::next(knot_first);
   }
   else if (time >= knot_last->first) {
       time = knot_last->first;
       prev = std::prev(knot_last);
       next = knot_last; 
   }
   else{ 
       next = knots.upper_bound(time);
       prev = std::prev(next);
   }
   
   // Assign p1 p2 values
   p1 = prev->second; t1 = prev->first;
   p2 = next->second; t2 = next->first;
  
   // "mirroring" strategy
   if (prev != knots.begin()) {
       p_prev = std::prev(prev);
       p0 = p_prev->second; t0 = p_prev->first;
   }
   else {
       p0 = p1 - (p2 - p1); t0 = t1 - (t2 - t1);
   }
   if (next != knot_last) {
       n_next = std::next(next);
       p3 = n_next->second; t3 = n_next->first;
   }
   else {
       p3 = p2 + (p2 - p1); t3 = t2 + (t2 - t1);
   }
   
   // normalized time
   double norm_factor = 1./ (t2 - t1);
   double norm_t  = (time - t1) * norm_factor;
   double norm_t0 = (t0   - t1) * norm_factor,
          norm_t1 = 0.,
	  norm_t2 = 1.,
	  norm_t3 = (t3   - t1) * norm_factor;
   // Tangents at p1 and p2
   m1 = (p2 - p0) / (norm_t2 - norm_t0);
   m2 = (p3 - p1) / (norm_t3 - norm_t1);
   
   return cubicSplineUnitInterval(p1, p2, m1, m2, norm_t, derivative);
}

// Removes the knot closest to the given time,
//    within the given tolerance..
// returns true iff a knot was removed.
template <class T>
inline bool Spline<T>::removeKnot(double time, double tolerance )
{
   // Empty maps have no knots.
   if( knots.size() < 1 )
   {
      return false;
   }

   // Look up the first element > or = to time.
   typename std::map<double, T>::iterator t2_iter = knots.lower_bound(time);
   typename std::map<double, T>::iterator t1_iter;
   t1_iter = t2_iter;
   t1_iter--;

   if( t2_iter == knots.end() )
   {
      t2_iter = t1_iter;
   }

   // Handle tolerance bounds,
   // because we are working with floating point numbers.
   double t1 = (*t1_iter).first;
   double t2 = (*t2_iter).first;

   double d1 = fabs(t1 - time);
   double d2 = fabs(t2 - time);


   if(d1 < tolerance && d1 < d2)
   {
      knots.erase(t1_iter);
      return true;
   }

   if(d2 < tolerance && d2 < d1)
   {
      knots.erase(t2_iter);
      return t2;
   }

   return false;
}

// Sets the value of the spline at a given time (i.e., knot),
// creating a new knot at this time if necessary.
template <class T>
inline void Spline<T>::setValue( double time, T value )
{
   knots[ time ] = value;
}

template <class T>
inline T Spline<T>::operator()( double time )
{
   return evaluate( time );
}
