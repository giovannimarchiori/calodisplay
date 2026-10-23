#include "PathUtils.h"

#include <cstdlib>
#include <filesystem>

namespace fs = std::filesystem;

namespace
{
std::string installDirectory;

// path of the running executable, empty if it cannot be determined
fs::path executablePath(const char* argv0)
{
  std::error_code ec;
  // Linux
  fs::path exe = fs::read_symlink("/proc/self/exe", ec);
  if (!ec)
    return exe;
  // otherwise argv[0], if it is a path
  if (argv0 && std::string(argv0).find('/') != std::string::npos)
  {
    exe = fs::canonical(argv0, ec);
    if (!ec)
      return exe;
  }
  return fs::path();
}
}

void PathUtils::setInstallDir(const char* argv0)
{
  const char* env = std::getenv("CALODISPLAY_DIR");
  if (env && *env)
  {
    installDirectory = env;
    return;
  }
  // the executable is in <installDir>/dist
  fs::path exe = executablePath(argv0);
  if (!exe.empty())
    installDirectory = exe.parent_path().parent_path().string();
  else
    installDirectory = fs::current_path().string();
}

const std::string& PathUtils::installDir()
{
  return installDirectory;
}

std::string PathUtils::findFile(const std::string& path)
{
  if (path.empty() || fs::exists(path) || fs::path(path).is_absolute() || installDirectory.empty())
    return path;
  fs::path inInstallDir = fs::path(installDirectory) / path;
  if (fs::exists(inInstallDir))
    return inInstallDir.string();
  return path;
}
