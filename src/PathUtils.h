#ifndef __PATHUTILS__
#define __PATHUTILS__

#include <string>

namespace PathUtils
{
  // set the directory where calodisplay is installed (the one containing
  // dist/, icons/, data/): $CALODISPLAY_DIR if set, otherwise the parent of
  // the directory of the executable
  void setInstallDir(const char* argv0);

  const std::string& installDir();

  // find a file distributed with calodisplay (icons, field map, geometry,
  // configuration): the path as given (absolute or relative to the current
  // directory) if it exists, otherwise relative to the installation directory.
  // Returns the path unchanged if it is found in neither
  std::string findFile(const std::string& path);
}

#endif
