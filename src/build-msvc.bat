set BuildDir=Build
set ExeName=dspsl.exe

set CFLAGS=/W3 /Ox /Ob2 /Oi /Ot /Oy /GL /D "WIN32" /D "NDEBUG" /D "_WINDOWS" /D "_UNICODE" /D "UNICODE" /FD /MT /GS- /Gy /arch:SSE2 /Zc:wchar_t- /GR- /Fo"%BuildDir%\\" /W3 /nologo /c /Zi /Gr /TC /Fd"%BuildDir%\\"

set LDFLAGS=/OUT:"%BuildDir%\%ExeName%" /NOLOGO /MANIFESTUAC:NO /IGNOREIDL /SUBSYSTEM:WINDOWS /OPT:REF /OPT:ICF /LTCG /RELEASE /DYNAMICBASE /NXCOMPAT:NO kernel32.lib user32.lib gdi32.lib comctl32.lib

mkdir %BuildDir%

rc /Fo "%BuildDir%\win32.res" ui\win32.rc
cl.exe %CFLAGS% mom.c mesh.c md.c global.c ui\win32.c
link.exe %LDFLAGS% build\*.obj build\*.res
