$ErrorActionPreference = "Stop"

$root = Split-Path -Parent (Split-Path -Parent $PSScriptRoot)
$buildDir = Join-Path $root ".bench_pybuild_variants"
$siteDir = Join-Path $root ".bench_resume_site"
$pybindDir = & "C:\Dev\Python\Python312\python.exe" -c "import pybind11; print(pybind11.get_cmake_dir())"

cmake -S (Join-Path $PSScriptRoot "native_variants") `
    -B $buildDir `
    -G "Visual Studio 17 2022" `
    -A x64 `
    -Dpybind11_DIR="$pybindDir" `
    -DPython_EXECUTABLE="C:\Dev\Python\Python312\python.exe"

cmake --build $buildDir --config Release

New-Item -ItemType Directory -Force $siteDir | Out-Null

Copy-Item (Join-Path $buildDir "Release\bench_variants_core.cp312-win_amd64.pyd") `
    (Join-Path $siteDir "bench_variants_core.cp312-win_amd64.pyd") -Force
