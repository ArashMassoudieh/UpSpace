#pragma once
#include <string>
#include <vector>
#include <map>

using namespace std;
class _command
{
public:
	_command();
	_command(const _command &c);
	_command &operator=(_command&);
	~_command();
	string command;
	map<string,string> parameters;
		
};

string insert_counter_in_file_name(string filename, int i);

