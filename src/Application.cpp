#include "Application.h"

#include <filesystem>
#include <stdexcept>
#include <thread>

constexpr static std::size_t MIN_PAR = 4;

constexpr static std::size_t POS_CFP = 1;
constexpr static std::size_t POS_UFP = 2;
constexpr static std::size_t POS_NFP = 3;
constexpr static std::size_t POS_OFP = 4;
constexpr static std::size_t POS_NOT = 5;
constexpr static std::size_t POS_WGR = 6;

Application::Application(gsl::span<char const *> const &args) : m_threadCount(std::thread::hardware_concurrency()) {
  if (args.size() >= MIN_PAR + 1) {
    readParameters(args);
  } else {
    throw std::runtime_error("Invalid parameter count");
  }
}

auto Application::checkIntegrity() const -> bool {
  std::filesystem::path contigsPath = m_contigsFilePath;
  std::filesystem::path unitigsPath = m_unitigsFilePath;
  std::filesystem::path nanoporePath = m_nanoporeFilePath;
  std::filesystem::path outputPath = m_outputPath;

  auto const exists = std::filesystem::exists(contigsPath) && std::filesystem::exists(unitigsPath) &&
                      std::filesystem::exists(nanoporePath) && std::filesystem::exists(outputPath);

  return exists;
}

void Application::readParameters(gsl::span<char const *> const &args) {

  m_contigsFilePath = args[POS_CFP];
  m_unitigsFilePath = args[POS_UFP];
  m_nanoporeFilePath = args[POS_NFP];
  m_outputPath = args[POS_OFP];

  // Check for optional parameter
  auto argIdx = args.size() > POS_WGR ? POS_WGR : args.size() - 1;
  switch (argIdx) {
  case POS_WGR:
    m_wiggleRoom = std::stoi(args[POS_WGR]);
    [[fallthrough]];
  case POS_NOT:
    m_threadCount = std::stoi(args[POS_NOT]);
    break;
  }
}
