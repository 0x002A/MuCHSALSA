#include "Application.h"

#include <filesystem>
#include <stdexcept>
#include <thread>

constexpr std::size_t MIN_PAR = 4;

constexpr std::size_t POS_CFP = 1;
constexpr std::size_t POS_UFP = 2;
constexpr std::size_t POS_NFP = 3;
constexpr std::size_t POS_OFP = 4;
constexpr std::size_t POS_NOT = 5;

Application::Application(const gsl::span<char *> &args) : m_threadCount(std::thread::hardware_concurrency()) {
  if (args.size() >= MIN_PAR + 1) {
    readParameters(args);
  } else {
    throw std::runtime_error("Invalid parameter count");
  }
}

bool Application::checkIntegrity() const {
  std::filesystem::path contigsPath = m_contigsFilePath;
  std::filesystem::path unitigsPath = m_unitigsFilePath;
  std::filesystem::path nanoporePath = m_nanoporeFilePath;
  std::filesystem::path outputPath = m_outputPath;

  const auto exists = std::filesystem::exists(contigsPath) && std::filesystem::exists(unitigsPath) &&
                      std::filesystem::exists(nanoporePath) && std::filesystem::exists(outputPath);

  return exists;
}

void Application::readParameters(const gsl::span<char *> &args) {

  m_contigsFilePath = args[POS_CFP];
  m_unitigsFilePath = args[POS_UFP];
  m_nanoporeFilePath = args[POS_NFP];
  m_outputPath = args[POS_OFP];

  // Check for optional parameter
  if (args.size() == MIN_PAR + 2) {
    m_threadCount = std::stoi(args[POS_NOT]);
  }
}
