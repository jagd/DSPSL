#include "global.h"

static void non_trace(TCHAR *s) {}

void (*mom_trace)(TCHAR *) = non_trace;
void (*mom_error)(TCHAR *) = non_trace;
