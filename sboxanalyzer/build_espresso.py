"""
Build espresso from source using CMake.
"""

import os
import shutil
import subprocess
import sys
from pathlib import Path


def build():
    """
    Build espresso from source code.
    
    Returns:
        str: Path to built espresso executable
        
    Raises:
        RuntimeError: If build fails
    """
    package_dir = Path(__file__).parent
    espresso_src = package_dir / "espresso_src"
    build_dir = package_dir / "espresso_build"
    
    if not espresso_src.exists():
        raise RuntimeError(f"Espresso source not found at {espresso_src}")
    
    # Check for CMake
    if not shutil.which("cmake"):
        raise RuntimeError(
            "CMake is required to build espresso.\n"
            "Install it with:\n"
            "  - macOS: brew install cmake\n"
            "  - Linux: sudo apt install cmake\n"
            "  - Windows: https://cmake.org/download/"
        )
    
    # Create build directory
    build_dir.mkdir(exist_ok=True)
    
    print(f"Building espresso in {build_dir}...")
    
    # Run CMake configure
    try:
        subprocess.run(
            ["cmake", "-B", str(build_dir), "-S", str(espresso_src)],
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"CMake configure failed:\n{e.stderr}")
    
    # Run CMake build
    try:
        subprocess.run(
            ["cmake", "--build", str(build_dir), "--config", "Release"],
            check=True,
            capture_output=True,
            text=True
        )
    except subprocess.CalledProcessError as e:
        raise RuntimeError(f"CMake build failed:\n{e.stderr}")
    
    # Find the built executable
    if sys.platform == "win32":
        espresso_exe = build_dir / "Release" / "espresso.exe"
        if not espresso_exe.exists():
            espresso_exe = build_dir / "espresso.exe"
    else:
        espresso_exe = build_dir / "espresso"
    
    if not espresso_exe.exists():
        raise RuntimeError(f"Built executable not found at {espresso_exe}")
    
    # Make executable on Unix
    if sys.platform != "win32":
        os.chmod(espresso_exe, 0o755)
    
    print(f"Successfully built espresso at {espresso_exe}")
    return str(espresso_exe)


if __name__ == "__main__":
    """Allow manual building with: python -m sboxanalyzer.build_espresso"""
    try:
        exe_path = build()
        print(f"\n✓ Espresso built successfully: {exe_path}")
    except Exception as e:
        print(f"\n✗ Build failed: {e}", file=sys.stderr)
        sys.exit(1)
