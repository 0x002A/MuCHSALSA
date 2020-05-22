#include <thread>
#include <stdexcept>
#include <filesystem>

#include "Application.h"

Application::Application(int argc, char* argv[])
  : m_contigsFilePath(nullptr)
  , m_unitigsFilePath(nullptr)
  , m_nanoporeFilePath(nullptr)
  , m_outputPath(nullptr)
  , m_threadCount(std::thread::hardware_concurrency())
{
  if(argc >= 4) {
    readParameters(argc, argv);
  } else {
    throw std::runtime_error("Invalid parameter count");
  }
}

bool
Application::checkIntegrity() const
{
  std::filesystem::path p_contigs = m_contigsFilePath;
  std::filesystem::path p_unitigs = m_unitigsFilePath;
  std::filesystem::path p_nanopore = m_nanoporeFilePath;
  std::filesystem::path p_output = m_outputPath;

  const auto exists = std::filesystem::exists(p_contigs)
                      && std::filesystem::exists(p_unitigs)
                      && std::filesystem::exists(p_nanopore)
                      && std::filesystem::exists(p_output);

  return exists;
}

void
Application::readParameters(int argc, char* argv[])
{
  m_contigsFilePath = argv[1];
  m_unitigsFilePath = argv[2];
  m_nanoporeFilePath = argv[3];
  m_outputPath = argv[4];

  // Optional one
  if (argc == 6) {
    m_threadCount = std::atoi(argv[5]);
  }
}
