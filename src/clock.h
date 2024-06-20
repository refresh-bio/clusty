/** @file */

#ifndef CLOCK_H
#define CLOCK_H


#include <chrono>

namespace ksi
{
   /** This class simplifies measurement of time. 
     @date 2018-04-07
     */
   class clock
   {
      std::chrono::high_resolution_clock _zegar;
      std::chrono::high_resolution_clock::time_point _start {};
      std::chrono::high_resolution_clock::time_point _stop  {};

      public:

      /** Starts the clock.
       * @return start time point
       @date 2018-04-07 */
      std::chrono::high_resolution_clock::time_point start ()
      {
         return _start = _zegar.now();
      }

      /** Stops the clock.
       * @return stop time point 
       @date 2018-04-07 
       */
      std::chrono::high_resolution_clock::time_point stop ()
      {
         return _stop = _zegar.now();
      }


      /** @return elapsed seconds 
       * @date 2018-04-07 
       */
      std::size_t elapsed_seconds()
      {
         return (std::chrono::duration_cast<std::chrono::seconds>(_stop - _start)).count();
      }


      /** @return elapsed milliseconds 
       * @date 2019-03-10 
       */
      std::size_t elapsed_milliseconds()
      {
         return (std::chrono::duration_cast<std::chrono::milliseconds>(_stop - _start)).count();
      }


   };
}

#endif 
