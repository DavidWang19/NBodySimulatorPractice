#include "OutputUtils.h"

#include <iomanip>

using namespace nbs;

OutputUtils::OutputUtils(const std::string& file_name) {
    m_output.open(file_name);
    m_output.setf(std::ios::fixed, std::ios::floatfield);
    m_output.precision(8);
}

OutputUtils::~OutputUtils() {
    m_output.close();
}

void Logger::log(std::string_view message) {
    auto now = std::chrono::system_clock::now();
    auto duration = std::chrono::duration_cast<std::chrono::milliseconds>(now - m_last_time);
    m_last_time = now;
    auto time = std::chrono::system_clock::to_time_t(now);
    auto time_str = std::put_time(std::localtime(&time), "%Y-%m-%d %H:%M:%S");
    m_output << "[" << time_str << "] " << message << " (" << duration.count() << "ms)\n";
}