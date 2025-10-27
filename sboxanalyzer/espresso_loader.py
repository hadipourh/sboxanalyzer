"""
Smart loader for espresso binary with automatic fallback to source build.
"""

import os
import sys
import platform
import subprocess
from pathlib import Path


def get_espresso_path():
    """
    Get the path to the espresso executable.
    
    Tries pre-built binary first, falls back to building from source if needed.
    
    Returns:
        str: Path to espresso executable
        
    Raises:
        RuntimeError: If espresso cannot be found or built
    """
    # Try pre-built binary first
    try:
        binary_path = _get_prebuilt_binary()
        if binary_path and _test_binary(binary_path):
            return binary_path
    except Exception as e:
        print(f"Pre-built binary not available: {e}")
    
    # Fallback: build from source
    print("Building espresso from source (this will take ~30 seconds, only once)...")
    try:
        built_path = _build_from_source()
        if built_path and _test_binary(built_path):
            return built_path
    except Exception as e:
        raise RuntimeError(
            f"Failed to build espresso from source: {e}\n"
            "Please ensure CMake and a C compiler are installed.\n"
            "See: https://github.com/hadipourh/espresso for build instructions"
        )
    
    raise RuntimeError("Could not find or build espresso executable")


def _get_prebuilt_binary():
    """Find and return path to pre-built binary for current platform."""
    package_dir = Path(__file__).parent
    binaries_dir = package_dir / "binaries"
    
    # Detect platform
    system = platform.system().lower()
    machine = platform.machine().lower()
    
    # Map to binary filename
    binary_map = {
        ("darwin", "arm64"): "espresso-darwin-arm64",
        ("darwin", "x86_64"): "espresso-darwin-x86_64",
        ("linux", "x86_64"): "espresso-linux-x86_64",
        ("linux", "aarch64"): "espresso-linux-aarch64",
        ("windows", "amd64"): "espresso-windows-x86_64.exe",
        ("windows", "x86_64"): "espresso-windows-x86_64.exe",
    }
    
    binary_name = binary_map.get((system, machine))
    if not binary_name:
        raise RuntimeError(f"No pre-built binary for {system}/{machine}")
    
    binary_path = binaries_dir / binary_name
    if not binary_path.exists():
        raise RuntimeError(f"Binary not found: {binary_path}")
    
    # Make executable on Unix
    if system != "windows":
        os.chmod(binary_path, 0o755)
    
    return str(binary_path)


def _test_binary(binary_path):
    """Test if binary works by running it with --help."""
    try:
        result = subprocess.run(
            [binary_path],
            input=b".i 0\n.o 0\n.e\n",
            capture_output=True,
            timeout=5
        )
        return result.returncode == 0
    except Exception:
        return False


def _build_from_source():
    """Build espresso from source using CMake."""
    from . import build_espresso
    return build_espresso.build()


# Cache the espresso path
_ESPRESSO_PATH = None

def get_espresso_path_cached():
    """Get espresso path with caching to avoid repeated checks."""
    global _ESPRESSO_PATH
    if _ESPRESSO_PATH is None:
        _ESPRESSO_PATH = get_espresso_path()
    return _ESPRESSO_PATH
