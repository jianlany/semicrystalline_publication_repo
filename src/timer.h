#pragma once
#include "string_tools.h"
#include <time.h>


class Timer {
    //static const int MY_CLOCK = CLOCK_THREAD_CPUTIME_ID;
    static const int MY_CLOCK = CLOCK_REALTIME;
public:
    Timer() {
        clock_gettime(MY_CLOCK, &_t0);
    }
    double elapsed() const {
        timespec _t1;
        clock_gettime(MY_CLOCK, &_t1);
        double dt = double(_t1.tv_sec - _t0.tv_sec);
        dt += double(_t1.tv_nsec - _t0.tv_nsec) * 1e-9;
        return dt;
    }
    std::string elapsed_str() const {
        double seconds = elapsed();
        if(seconds < 0.0) return make_string(seconds)+" ss";
        unsigned t = unsigned(seconds);
        string time, time_fmt;
        if(t>=3600){
            time += make_string(t/3600,2) + ":";
            time_fmt += "hh:";
            t %= 3600;
        }
        if(t>=60 || !time_fmt.empty()){
            time += make_string(t/60,2) + ":";
            time_fmt += "mm:";
            t %= 60;
        }
        // preserve decimal point
        seconds -= floor(seconds/60.0)*60.0;
        if(seconds<10.0) time += "0";
        time += make_string(seconds, 4, 2);
        time_fmt += "ss";
        return time + " " + time_fmt;
    }

private:
    timespec _t0;
};

