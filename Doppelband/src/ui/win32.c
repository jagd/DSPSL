#include <windows.h>
#include <stdio.h>
#include <math.h>

#include "resource.h"

#include "../global.h"
#include "../mom.h"

#define VERSION_STR TEXT("0.3")

#define ABOUT_INFO \
	TEXT("Calculate the characteristic Impedance ")\
	TEXT("of a double strip\ntransmission line ")\
	TEXT("using the Method of Moments.\n\n")\
	TEXT("All of the calculation libraries ")\
	TEXT("(including the matrix operation)\n")\
	TEXT("are wrote by myself. ")\
	TEXT("The codes are of course for this problem\n")\
	TEXT("heavily optimized.\n\n")\
	TEXT("Email : chuanren.wu@student.kit.edu")


#define TEXT_BUF_LENGTH 1000
TCHAR buf[TEXT_BUF_LENGTH + 1];
double w1, w2, d, h, eps_r, port_ext;

HWND hMainDlg;
HINSTANCE hInst;

/* n > 0*/
static double nocrt_pow(double x,int n)
{
	double result;

	if (n == 0) {
		return 1;
	} else if (n < 0) {
		return nocrt_pow(1/x, -n);
	} else if (n == 1) {
		return x;
	}

	result = nocrt_pow(x, (n>>1));
	result *= result;

	if(n & 1){ /*wenn nichtgerade Zahl*/
		result*=x;
	}

	return result;
}

static void Chop(TCHAR* str)
{
	int i, j;

	/* chop the head */
	i = 0;
	while (str[i] == ' '
		|| str[i] == '\t'
		|| str[i] == '\n'
		|| str[i] == '\r'
		) {
 		++i;
	}

	j = 0;
	while (str[i] != 0) {
		str[j++] = str[i++];
	}
	str[j] = 0;

	/* chop the tail */
	--j;
	while (j >= 0
		&& ( str[j] == ' '
		  || str[j] == '\t'
		  || str[j] == '\n'
		  || str[j] == '\r'
		  )) {
		--j;
	}
	str[j+1] = 0;
}

/*
   (sw)scanf do not support formal language check.
   So implement myself
   The string will be in-place chopped
*/
static int StringToDouble(TCHAR* str, double *val)
{
	int i = 0;
	int sign = 1, esign = 1;
	double x = 0;
	int xe = 0, e = 0;
	char flag = 1;

	Chop(str);

	if (str[i] == (TCHAR)'-') {
		sign = -1;
		++i;
	}

	while (str[i] >= (TCHAR)'0' && str[i] <= (TCHAR)'9') {
		x *= 10;
		x += str[i++] - (TCHAR)'0';
		flag = 0;
	}

	if (sign < 0) {
		x = -x;
	}

	if (str[i] == (TCHAR)'.') {
		++i;
		while (str[i] >= (TCHAR)'0' && str[i] <= (TCHAR)'9') {
			x *= 10;
			x += str[i++] - (TCHAR)'0';
			++xe;
		}
		flag = 0;
	}

	if (flag) {
		return -1;
	}

	if (str[i] == 'e' || str[i] == 'E') {
		flag = 1;

		++i;
		if (str[i] == (TCHAR)'-') {
			esign = -1;
			++i;
		}
		while (str[i] >= (TCHAR)'0' && str[i] <= (TCHAR)'9') {
			e *= 10;
			e += str[i++] - (TCHAR)'0';
			flag = 0;
		}

		if (esign < 0) {
			e = -e;
		}
	}

	if (flag) {
		return -1;
	}

	e -= xe;

	if (str[i] != 0) {
		return -1; /* error */
	} else {
		*val = x * nocrt_pow(10, e);
		return 0;
	}
}

static void trace(LPTSTR msg)
{
	int n;
	SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_ADDSTRING, 0, (LPARAM)msg);
	n = SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_GETCOUNT, 0, 0);
	SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_SETCURSEL, n-1, 0);
}


DWORD WINAPI Calc(LPVOID lpParam)
{
	struct MeshConfig *conf;
	struct MD *x[2], *x0[2];
	double c[2];
	double z0;

	EnableWindow(lpParam, FALSE);

	/* calc */
	conf = mesh_new(
		w1, w2,
		d,
		port_ext,
		h,
		eps_r
		);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("Mesh cells = %d"), conf->index[ID_MESH_CELLS]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("\tcells for strips = %d"), conf->index[ID_STRIP_END]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("\tcells for dielectrics = %d"),
			conf->index[ID_DIELECTRIC_END]
				-conf->index[ID_DIELECTRIC_START]);
	trace(buf);

	z0 = mom(conf, x0, x, c);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("Effective permittivity = %lf"), c[1]/c[0]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH, TEXT("Z0 = %lf Ohm"), z0);
	trace(buf);

	md_free(x[0]);
	md_free(x[1]);
	md_free(x0[0]);
	md_free(x0[1]);

	mesh_free(conf);

	EnableWindow(lpParam, TRUE);

	return 0;
}


int Read(HWND hDlg)
{
	SendDlgItemMessage(hDlg, IDC_LIST_RESULT, LB_RESETCONTENT, 0, 0);
	SendDlgItemMessage(hDlg, IDC_LIST_RESULT, LB_ADDSTRING, 0,
		(LPARAM)TEXT("\tCalculator version ")VERSION_STR);


	GetDlgItemText(hDlg, IDC_EDIT_W1, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &w1)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}

	w1 *= 1e-3;
	if (w1 < 1e-20) {
		mom_error(TEXT("ERROR: The value of `w1` must be > 0"));
		return 1;
	}


	GetDlgItemText(hDlg, IDC_EDIT_W2, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &w2)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}
	if (w2 < 1e-20) {
		mom_error(TEXT("ERROR: The value of `w2` must be > 0"));
		return 1;
	}
	w2 *= 1e-3;


	GetDlgItemText(hDlg, IDC_EDIT_H, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &h)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}
	if (h < 1e-20) {
		mom_error(TEXT("ERROR: The value of `h` must be > 0"));
		return 1;
	}
	h *= 1e-3;


	GetDlgItemText(hDlg, IDC_EDIT_PORT, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &port_ext)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}
	if (port_ext < 1e-20) {
		mom_error(TEXT("small `p` will be supported in the next version"));
		return 1;
	}
	port_ext *= 1e-3;


	GetDlgItemText(hDlg, IDC_EDIT_EPS, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &eps_r)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}
	if (eps_r < 1.0) {
		mom_error(TEXT("ERROR: epsilon_r must be > 1"));
		return 1;
	}


	GetDlgItemText(hDlg, IDC_EDIT_D, buf, TEXT_BUF_LENGTH);
	if (0 != StringToDouble(buf, &d)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	}
	d *= 1e-3;


	return 0;
}


INT_PTR CALLBACK MainWndProc(
	HWND hDlg,
	UINT uMsg,
	WPARAM wParam, LPARAM lParam)
{
	DWORD dwTID;

	switch (uMsg) {
		case WM_INITDIALOG:
			hMainDlg = hDlg;
			SetDlgItemText(hDlg, IDC_EDIT_W1, TEXT("2.4"));
			SetDlgItemText(hDlg, IDC_EDIT_W2, TEXT("10"));
			SetDlgItemText(hDlg, IDC_EDIT_H, TEXT("0.79"));
			SetDlgItemText(hDlg, IDC_EDIT_PORT, TEXT("2.0"));
			SetDlgItemText(hDlg, IDC_EDIT_EPS, TEXT("2.2"));
			SetDlgItemText(hDlg, IDC_EDIT_D, TEXT("0.0"));
			SetDlgItemText(hDlg, IDC_EDIT_MESHSTEP, TEXT("auto"));

			SetClassLong(hDlg, GCL_HICON, 
			    (LONG)LoadIcon(hInst, MAKEINTRESOURCE(IDI_ICON))); 
			break;
		case WM_COMMAND:
			switch (LOWORD(wParam)) {
				case IDOK:
					if (Read(hDlg) == 0) {
						CreateThread(
							NULL,
							0,
							Calc,
							GetDlgItem(hDlg, IDOK),
							0,
							&dwTID);
					}
					break;
				case IDHELP:
					MessageBox(
						hDlg,
						ABOUT_INFO,
						TEXT("ABOUT"),
						MB_OK);
				break;
			}
			break;
		case WM_CLOSE:
			EndDialog(hDlg, 0);
			break;
		default:
			return FALSE;
	}
	return TRUE;
}


static void inter_trace(TCHAR *msg)
{
	/* for non-unicode */
	/*
	OemToCharBuff(msg, buf, TEXT_BUF_LENGTH);
	trace(buf);
	 */

	trace(msg);
}

int WINAPI WinMain(HINSTANCE hInstance,
		   HINSTANCE hPrevInstance,
		   LPSTR lpCmdLine,
		   int nShowCmd)
{
	hInst = hInstance;

	mom_trace = inter_trace;
	mom_error = mom_trace;

	return DialogBoxParam(
			hInstance,
			MAKEINTRESOURCE(IDD_DLG_MAIN),
			HWND_DESKTOP,
			MainWndProc,
			0
		);
}
