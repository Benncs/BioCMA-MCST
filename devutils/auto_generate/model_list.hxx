#ifndef __GENERATED_MODEL_LIST__HPP__
#define __GENERATED_MODEL_LIST__HPP__

#include <models/types.hpp>
#include <vector> 

KModel load_model_(const std::string& name,bool use_python_module=true);


std::vector<std::string> get_available_models(bool use_python_module=true);

#endif 
