#include "Application.h"

#include <filesystem>
#include <gsl/span>
#include <stdexcept>
#include <thread>

Application::Application(int argc, char *argv[]) : m_threadCount(std::thread::hardware_concurrency()) {
  if (argc >= 4) {
    readParameters(argc, argv);
  } else {
    throw std::runtime_error("Invalid parameter count");
  }
}

bool Application::checkIntegrity() const {
  std::filesystem::path p_contigs = m_contigsFilePath;
  std::filesystem::path p_unitigs = m_unitigsFilePath;
  std::filesystem::path p_nanopore = m_nanoporeFilePath;
  std::filesystem::path p_output = m_outputPath;

  const auto exists = std::filesystem::exists(p_contigs) && std::filesystem::exists(p_unitigs) &&
                      std::filesystem::exists(p_nanopore) && std::filesystem::exists(p_output);

  return exists;
}

void Application::readParameters(int argc, char *argv[]) {
  gsl::span<char *> args = {argv, static_cast<std::size_t>(argc)};

  m_contigsFilePath = args[1];
  m_unitigsFilePath = args[2];
  m_nanoporeFilePath = args[3];
  m_outputPath = args[4];

  // Optional one
  if (argc == 6) {
    m_threadCount = std::stoi(args[5]);
  }
}
