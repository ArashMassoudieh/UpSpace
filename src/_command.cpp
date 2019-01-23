#include "_command.h"
#include "StringOP.h"



_command::_command()
{
	parameters["z_factor"] = "1";
	parameters["offset"] = "0";
	parameters["nbins"] = "10";
	parameters["log"] = "0";
	parameters["color"] = "0";
	parameters["smoothing_factor"] = "0";

}

_command::_command(const _command & c)
{
	command = c.command;
	parameters = c.parameters;
}

_command & _command::operator=(_command &c)
{
	command = c.command;
	parameters = c.parameters;

	return *this;
}


_command::~_command()
{
}

string insert_counter_in_file_name(string filename, int i)
{
	string name = filename.substr(0, filename.size() - 4);
	string ext = filename.substr(filename.size() - 3, 3);
	string out = name + "_" + numbertostring(i) + "." + ext;
	return out;

}
