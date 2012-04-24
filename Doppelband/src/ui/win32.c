#include <windows.h>
#include <stdio.h>

#include "resource.h"

#include "../global.h"
#include "../mom_mesh.h"
#include "../mom.h"

#define VERSION_STR TEXT("0.2")


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


static void trace(LPTSTR msg)
{
	SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_ADDSTRING, 0, (LPARAM)msg);
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
		TEXT("Mesh cells = %d\n"), conf->index[ID_MESH_CELLS]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("\tcells for strips = %d\n"), conf->index[ID_STRIP_END]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("\tcells for dielectrics = %d\n"),
			conf->index[ID_DIELECTRIC_END]
				-conf->index[ID_DIELECTRIC_START]);
	trace(buf);

	z0 = mom(conf, x0, x, c);

	swprintf_s(buf, TEXT_BUF_LENGTH,
		TEXT("Effective permittivity = %lf\n"), c[1]/c[0]);
	trace(buf);

	swprintf_s(buf, TEXT_BUF_LENGTH, TEXT("Z0 = %lf Ohm\n"), z0);
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
	swscanf_s(buf, TEXT("%le"), &w1);
	w1 *= 1e-3;
	if (w1 < 1e-20) {
		mom_error(TEXT("ERROR: The value of `w1` must be > 0"));
		return 1;
	}

	GetDlgItemText(hDlg, IDC_EDIT_W2, buf, TEXT_BUF_LENGTH);
	swscanf_s(buf, TEXT("%le"), &w2);
	if (w2 < 1e-20) {
		mom_error(TEXT("ERROR: The value of `w2` must be > 0"));
		return 1;
	}
	w2 *= 1e-3;

	GetDlgItemText(hDlg, IDC_EDIT_H, buf, TEXT_BUF_LENGTH);
	swscanf_s(buf, TEXT("%le"), &h);
	if (h < 1e-20) {
		mom_error(TEXT("ERROR: The value of `h` must be > 0"));
		return 1;
	}
	h *= 1e-3;

	GetDlgItemText(hDlg, IDC_EDIT_PORT, buf, TEXT_BUF_LENGTH);
	swscanf_s(buf, TEXT("%le"), &port_ext);
	if (port_ext < 1e-20) {
		mom_error(TEXT("small `p` will be supported in the next version"));
		return 1;
	}
	port_ext *= 1e-3;

	GetDlgItemText(hDlg, IDC_EDIT_EPS, buf, TEXT_BUF_LENGTH);
	swscanf_s(buf, TEXT("%le"), &eps_r);
	if (eps_r < 1.0) {
		mom_error(TEXT("ERROR: epsilon_r must be > 1"));
		return 1;
	}

	GetDlgItemText(hDlg, IDC_EDIT_D, buf, TEXT_BUF_LENGTH);
	swscanf_s(buf, TEXT("%le"), &d);
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
