#include <iostream>
#include <stdexcept>
#include <string>

// mock for UT only
namespace ModuleBase
{
void WARNING(const std::string &file,const std::string &description)
{
    return;
}

void WARNING_QUIT(const std::string &file,const std::string &description, int ret)
{
#ifdef __NORMAL

		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;

#else
		std::cout << " " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << " " << std::endl;
		std::cout << " " << description << std::endl;
		std::cout << " " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
		std::cout << "                         NOTICE                           " << std::endl;
		std::cout << " !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!" << std::endl;
#endif

    throw std::runtime_error("WARNING_QUIT");
}

void WARNING_QUIT(const std::string &file,const std::string &description)
{
	WARNING_QUIT(file, description, 0);
}

}
