#ifndef DRIVER_H
#define DRIVER_H

class Driver
{
  public:
    Driver();
    ~Driver();

    /**
     * @brief An interface function.
     * This function calls "this->reading()", "this->atomic_world()" in order.
     */
    void init();

  private:
    /**
     * @brief Reading the parameters and split the MPI world.
     * This function read the parameter in "INPUT", "STRU" etc,
     * and split the MPI world into different groups.
     */
    void reading();

    /**
     * @brief An interface function.
     * This function calls "this->driver_run()" to do calculation,
	 * and log the time and  memory consumed during calculation.
     */
    void atomic_world();

    // the actual calculations
    void driver_run();
};

#endif
