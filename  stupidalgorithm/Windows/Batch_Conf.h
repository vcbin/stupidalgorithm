#ifndef STUPIDALGO_BATCH_CONF_CLASS
#define STUPIDALGO_BATCH_CONF_CLASS

#include <string>
#include <vector>

struct task_info
{
	int algo_code;
	std::string conf_path;
};

typedef std::vector<task_info> vec_task;

struct batch_conf
{
public:
	// ctor
	batch_conf(std::string path=""):m_strpath(path){}

	int load_conf(std::string path="");
	std::string m_strpath;
	vec_task vec_task;
};

#endif