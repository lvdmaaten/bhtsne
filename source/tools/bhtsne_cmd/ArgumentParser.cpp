#include "ArgumentParser.h"

#include <iostream>

namespace cppassist
{
    ArgumentParser::ArgumentParser() = default;

    ArgumentParser::~ArgumentParser() = default;

    void ArgumentParser::parse(int argc, char * argv[])
    {
        m_options.clear();
        m_params.clear();
        auto argumentVector = std::vector<std::string>(argv, argv + argc);

        for (auto it = argumentVector.begin() + 1; it != argumentVector.end(); ++it)
        {
            // Get current and next argument
            auto arg = *it;

            if (arg.find("--") == 0) // Options with value (--option-name <value>)
            {
                // Save value
                ++it;
                if (it != argumentVector.end()) {
                    m_options[arg] = *it;
                }
                else
                {
                    m_options[arg] = "";
                }
            }
            else if (arg.find('-') == 0) // Options without value (-option-name)
            {
                m_options[arg] = "true";
            }
            else // Additional parameters
            {
                m_params.push_back(arg);
            }
        }
    }

    const std::map<std::string, std::string> & ArgumentParser::options() const
    {
        return m_options;
    }

    bool ArgumentParser::isSet(const std::string & option) const
    {
        return m_options.count(option) > 0;
    }

    const std::vector<std::string> & ArgumentParser::params() const
    {
        return m_params;
    }
} // namespace cppassist
