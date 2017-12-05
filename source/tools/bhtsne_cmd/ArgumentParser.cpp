#include "ArgumentParser.h"

#include <iostream>


namespace cppassist
{


    ArgumentParser::ArgumentParser()
    {
    }

    ArgumentParser::~ArgumentParser()
    {
    }

    void ArgumentParser::parse(int argc, char * argv[])
    {
        m_options.clear();
        m_params.clear();

        for (int i=1; i<argc; i++)
        {
            // Get current and next argument
            std::string arg  = argv[i];
            std::string next = (i+1 < argc ? argv[i+1] : "");

            // Options with value (--option-name <value>)
            if (arg.find("--") == 0)
            {
                // Save value
                m_options[std::move(arg)] = std::move(next);
                i++;
            }

                // Options without value (-option-name)
            else if (arg.find("-") == 0)
            {
                m_options[std::move(arg)] = "true";
            }

                // Additional parameters
            else
            {
                m_params.push_back(std::move(arg));
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