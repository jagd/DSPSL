set BuildDir=Build
set ExeName=win32.exe

set CFLAGS=-O3 -o%BuildDir%\\%ExeName% -DUNICODE -DMOM_ENABLE_COMCTL -D_WIN32_IE=0x0500 -D_WIN32_WINNT=0x0501 -Dswprintf_s=snwprintf

mkdir %BuildDir%

windres -o "%BuildDir%\resource.o" ui\win32.rc
windres -o "%BuildDir%\manifest.o" ui\win32_manifest.rc

gcc %CFLAGS% -mwindows mom.c mesh.c md.c global.c ui\win32.c %BuildDir%\resource.o %BuildDir%\manifest.o -lcomctl32 -lgdi32
strip %BuildDir%\\%ExeName%
