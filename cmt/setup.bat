@echo off
if .%1==. (set tag=VC7debug ) else set tag=%1
set tempfile=%HOME%\tmpsetup.bat
c:\Programs\CMT\v1r14p20031120\VisualC\cmt.exe -quiet -bat -pack=calibGenTKR -version=v0r0 setup -tag=%tag% >"%tempfile%"
if exist "%tempfile%" call "%tempfile%"
if exist "%tempfile%" del "%tempfile%"
set PATH=%LD_LIBRARY_PATH%;%PATH%