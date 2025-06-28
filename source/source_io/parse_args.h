#ifndef ModuleIO_PARSE_ARGS_H
#define ModuleIO_PARSE_ARGS_H

namespace ModuleIO
{
/**
 * @brief This function reture the version information when using command
 * "abacus --version", "abacus -v" or "abacus -V". Otherwise, it does nothing.
 * 
 * @param [in] argc (ARGument Count) is an integer variable that stores the number 
 * of command-line arguments passed by the user including the name of the program. 
 * So if we pass a value to a program, the value of argc would be 2 (one for 
 * argument and one for program name) 
 * @param [in] argv (ARGument Vector) is an array of character pointers listing
 * all the arguments. If argc is greater than zero, the array elements from 
 * argv[0] to argv[argc-1] will contain pointers to strings. argv[0] is the name 
 * of the program , After that till argv[argc-1] every element is command -line 
 * arguments.
 */
void parse_args(int argc, char** argv);
} // namespace ModuleIO

#endif