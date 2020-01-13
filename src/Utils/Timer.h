/* 
* Copyright (C) 2015, Adrian Jarabo (http://giga.cps.unizar.es/~ajarabo/)
*
* Permission is hereby granted, free of charge, to any person obtaining
* a copy of this software and associated documentation files (the "Software"),
* to deal in the Software without restriction, including without limitation
* the rights to use, copy, modify, merge, publish, distribute, sublicense,
* and/or sell copies of the Software, and to permit persons to whom
* the Software is furnished to do so, subject to the following conditions:

* The above copyright notice and this permission notice shall be included
* in all copies or substantial portions of the Software.

* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
* EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF
* MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.
* IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM,
* DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT,
* TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE
* OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.
*/

#ifndef __OSTIMER_H__
#define __OSTIMER_H__

// Created by bdl 5. april 2002
// The purpose of this file is to make a timer function that is as 
// precise as posible on any given platform

#include "bunnykiller.h"

#ifdef WIN32
#include <windows.h>
#include <time.h>
static LARGE_INTEGER	largeInteger;
#else
#include <sys/time.h>
#endif

//namespace Util
//{
		/** A simple timer for measuring "wall clock" time - i.e. time 
				as it passes on a wall clock and not, say, clock ticks consumed by
				the program. */
	class Timer
	{
#ifdef WIN32
		double freq;
		double start_count;
#else
		timeval start_time;
#endif
	public:

		/// Start timer
		void start()
		{
#ifdef WIN32
			QueryPerformanceFrequency(&largeInteger);
			freq = static_cast<double>(largeInteger.QuadPart);
			QueryPerformanceCounter(&largeInteger);
			start_count = static_cast<double>(largeInteger.QuadPart);
#else
			gettimeofday(&start_time, nullptr);
#endif
		}
		
		/// Return number of seconds since start was called.
		double get_secs()
		{
#ifdef WIN32
			QueryPerformanceCounter(&largeInteger);
			double now_count = static_cast<double>(largeInteger.QuadPart);
			return (now_count-start_count)/freq;
#else
			timeval now;
			gettimeofday(&now, nullptr);
			return (now.tv_sec-start_time.tv_sec) + 
				(now.tv_usec-start_time.tv_usec)/1.0e6;
#endif
		}
	};
//}
#endif
