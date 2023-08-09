#ifndef TIMER_H
#define TIMER_H

#include <ctime>
#include <vector>

class Timer
{

private:
	std::clock_t c_ini;
	std::clock_t c_end;
	bool running = false;
	double time = 0.0;

public:
	void start(void) 
	{ 
		running = true;
		c_ini = std::clock(); 
	}
	void stop(void)
	{
		if (running)
		{
			c_end = std::clock();
			time = 1000.0*(c_end-c_ini)/CLOCKS_PER_SEC;
			running = false;
		}
	}
	double get_time(void) { return time; }

};

#endif /* TIMER_H */
