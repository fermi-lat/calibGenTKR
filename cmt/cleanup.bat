set tempfile=%HOMEDRIVE%%HOMEPATH%tmpsetup.bat
%CMTROOT%\%CMTBIN%\cmt.exe -quiet cleanup -bat -path=C:\Glast\em\BadStripsCalib %1 %2 %3 %4 %5 %6 %7 %8 %9 >%tempfile%
if exist %tempfile% call %tempfile%
if exist %tempfile% del %tempfile%

