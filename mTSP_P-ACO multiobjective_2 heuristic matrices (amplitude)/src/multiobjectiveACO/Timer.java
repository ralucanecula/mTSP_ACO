package multiobjectiveACO;

public class Timer {

    private static long startTime;

    //virtual and real time of day are computed and stored to allow at later time the computation of the elapsed time
    static void start_timers()
    {
    	startTime = System.currentTimeMillis();
    }

    //return the time used in seconds (virtual or real, depending on type)
    static double elapsed_time()
    {
    	return (System.currentTimeMillis() - startTime) / 1000.0;
    }

}
