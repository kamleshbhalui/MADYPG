#include "framesio.h"

#include <bitsery/adapter/stream.h>

#include <filesystem>
#include <fstream>

#include "../utils/debug_logging.h"

namespace fs = std::filesystem;

bool deserialize_frames(const std::string& filepath,
                        std::vector<Frame>& frames) {
  frames.clear();

  if (!fs::is_regular_file(filepath)) {
    Debug::error("Could not load binary file (not a file):", filepath);
    return false;
  }

  std::fstream fstrm(filepath, std::ios::binary | std::ios::in);
  if (!fstrm.is_open()) {
    Debug::error("Could not load binary file (fstream):", filepath);
    return false;
  }

  // first = error code, second = is buffer was successfully read from begin to
  // the end.
  auto state =
      bitsery::quickDeserialization<bitsery::InputStreamAdapter>(fstrm, frames);
  if (!Debug::msgassert(
          "Could not deserialize " + filepath,
          state.first == bitsery::ReaderError::NoError && state.second))
    return false;
  return true;
}