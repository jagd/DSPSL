#include <stdio.h>
#include <math.h>
#include <windows.h>

/* without commctrl is also usable , only for visual styles */
/* #define MOM_ENABLE_COMCTL */

#ifdef MOM_ENABLE_COMCTL
#include <commctrl.h>
#endif

#include "resource.h"

#include "../global.h"
#include "../mom.h"

#define VERSION_STR TEXT("0.4")

#define TEXT_BUF_LENGTH 1000
TCHAR buf[TEXT_BUF_LENGTH + 1];
double w1, w2, d, h, eps_r, port_ext, mesh_step;

struct PlotData {
	struct MeshConfig *conf;
	struct MD *x;
};

TCHAR StrInf[] = TEXT("INF");
TCHAR StrAuto[] = TEXT("AUTO");

HWND hMainDlg;
HINSTANCE hInst;
HBITMAP hLayoutBMP;
BITMAP LayoutBMP;

DWORD TIDCalc;
HANDLE hThreadCalc;

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
	LRESULT n;
	SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_ADDSTRING, 0, (LPARAM)msg);
	n = SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_GETCOUNT, 0, 0);
	SendDlgItemMessage(hMainDlg, IDC_LIST_RESULT,
		LB_SETCURSEL, n-1, 0);
}


void PlotCharge(HWND hWnd, struct PlotData *pd)
{
	HDC hdc;
	PAINTSTRUCT paint;
	RECT rtCanvas;
	int i;
	double x_max, x_min, y_max, y_min;
	/* the ratio between screen pixel and the mesh coordinate */
	double y_scale, x_scale;
	/* the original point in client window coordnate: */
	double o_hori, o_vert;
	HBRUSH brush[3];

	brush[0] = CreateSolidBrush(RGB(255,255,255));
	brush[1] = CreateSolidBrush(RGB(255,0,0));
	brush[2] = CreateSolidBrush(RGB(0,0,255));

	/*
		-------->
		|  hori
		|
		|
		|  vert
		v
	*/

	GetClientRect(hWnd,&rtCanvas);

	hdc = BeginPaint(hWnd, &paint);

	x_max = pd->conf->mesh[0].centre;
	x_min = pd->conf->mesh[0].centre;
	y_max = pd->x->buf[0];
	y_min = pd->x->buf[0];
	for (i = 0; i < pd->conf->index[ID_MESH_CELLS]; ++i) {
		double centre = pd->conf->mesh[i].centre;
		double hw = pd->conf->mesh[i].hw;
		double y = pd->x->buf[i];

		if (x_max < centre + hw) {
			x_max = centre + hw;
		}
		if (x_min > centre - hw) {
			x_min = centre - hw;
		}

		if (y > y_max) {
			y_max = y;
		}
		if (y < y_min) {
			y_min = y;
		}
	}
	x_scale = (double)(rtCanvas.right - rtCanvas.left) / (x_max - x_min);
	o_hori = rtCanvas.right/2 + (x_min + (x_max-x_min)/2) * x_scale;
	y_scale = (double)(rtCanvas.bottom - rtCanvas.top) / (y_max - y_min);
	o_vert = y_max * y_scale;

	for (i = pd->conf->index[ID_STRIP0_START];
		i < pd->conf->index[ID_STRIP0_END]; ++i) {
		double centre, hw;
		RECT rt;

		centre = pd->conf->mesh[i].centre;
		hw = pd->conf->mesh[i].hw;

		rt.top = (LONG)(o_vert - pd->x->buf[i]*y_scale);
		rt.bottom = (LONG)(o_vert) + 1; /* because the FillRect function */
		rt.left = (LONG)(o_hori + (centre-hw)*x_scale);
		rt.right = (LONG)(o_hori + (centre+hw)*x_scale) + 1;
		FillRect(hdc, &rt, brush[1]);
	}

	for (i = pd->conf->index[ID_DIELECTRIC0_START];
		i < pd->conf->index[ID_DIELECTRIC0_START]; ++i) {
		double centre, hw;
		RECT rt;

		centre = pd->conf->mesh[i].centre;
		hw = pd->conf->mesh[i].hw;

		rt.top = (LONG)(o_vert - pd->x->buf[i]*y_scale);
		rt.bottom = (LONG)(o_vert) + 1; /* because the FillRect function */
		rt.left = (LONG)(o_hori + (centre-hw)*x_scale);
		rt.right = (LONG)(o_hori + (centre+hw)*x_scale) + 1;
		FillRect(hdc, &rt, brush[1]);
	}

	for (i = pd->conf->index[ID_STRIP1_START];
		i < pd->conf->index[ID_STRIP1_END]; ++i) {
		double centre, hw;
		RECT rt;

		centre = pd->conf->mesh[i].centre;
		hw = pd->conf->mesh[i].hw;

		rt.top = (LONG)(o_vert - pd->x->buf[i]*y_scale);
		rt.bottom = (LONG)(o_vert) + 1; /* because the FillRect function */
		rt.left = (LONG)(o_hori + (centre-hw)*x_scale);
		rt.right = (LONG)(o_hori + (centre+hw)*x_scale) + 1;
		FillRect(hdc, &rt, brush[2]);
	}

	for (i = pd->conf->index[ID_DIELECTRIC1_START];
		i < pd->conf->index[ID_DIELECTRIC1_START]; ++i) {
		double centre, hw;
		RECT rt;

		centre = pd->conf->mesh[i].centre;
		hw = pd->conf->mesh[i].hw;

		rt.top = (LONG)(o_vert - pd->x->buf[i]*y_scale);
		rt.bottom = (LONG)(o_vert) + 1; /* because the FillRect function */
		rt.left = (LONG)(o_hori + (centre-hw)*x_scale);
		rt.right = (LONG)(o_hori + (centre+hw)*x_scale) + 1;
		FillRect(hdc, &rt, brush[2]);
	}

	EndPaint(hWnd, &paint);

	DeleteObject(brush[0]);
	DeleteObject(brush[1]);
	DeleteObject(brush[2]);
}


INT_PTR CALLBACK PlotWndProc(
	HWND hDlg,
	UINT uMsg,
	WPARAM wParam, LPARAM lParam)
{
	static struct PlotData *pd;

	switch (uMsg) {
		case WM_SIZE:
			InvalidateRect(hDlg, NULL, TRUE);
			break;
		case WM_INITDIALOG:
			pd = (struct PlotData*)lParam;
			break;
		case WM_PAINT:
			/* pd != NULL */
			if (GetUpdateRect(hDlg, NULL, FALSE) == 0) {
				break;
			}
			PlotCharge(hDlg, pd);
			break;
		case WM_CLOSE:
			EndDialog(hDlg, 0);
			break;
		default:
			return FALSE;
	}
	return TRUE;
}


void ModalPlotWindow(struct PlotData *pd)
{
	DialogBoxParam(hInst, MAKEINTRESOURCE(IDD_PLOT),
		hMainDlg, PlotWndProc, (LPARAM)pd);
}

DWORD WINAPI Calc(LPVOID lpParam)
{
	struct MeshConfig *conf;
	struct MD *x[2], *x0[2];
	double c[2];
	double z0;
	char flag_inf = 0;

	EnableWindow(GetDlgItem(hMainDlg, IDC_CALC), FALSE);

	/* calc */

	/* infinite w2 (fall back to microstrip) */
	if (w2 <= 0) {
		flag_inf = 1;
		w2 = w1;
		d = 0;
		h *= 2;

		trace(TEXT("Infinite w2 (fall back to microstrip)"));
		trace(TEXT("Calculation with two symmetric strips"));
	}

	/* auto port_ext */
	if (port_ext <= 0) {
		port_ext = 3.0*h;

		swprintf_s(buf, TEXT_BUF_LENGTH,
			TEXT("use p = %lf mm ")
			TEXT("(emulation of a infinite large p)"),
			1e3*port_ext);
		trace(buf);
	}

	conf = mesh_new(
		w1, w2,
		d,
		port_ext,
		h,
		eps_r,
		mesh_step
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

	swprintf_s(buf, TEXT_BUF_LENGTH, TEXT("Z0 = %lf Ohm"),
		flag_inf ? z0/2 : z0);
	trace(buf);

	if (BST_CHECKED == SendDlgItemMessage(hMainDlg,
		IDC_CHECK_PLOT, BM_GETSTATE, 0, 0)) {
			struct PlotData pd;
			pd.conf = conf;
			pd.x = x[0];
			ModalPlotWindow(&pd);
	}

	md_free(x[0]);
	md_free(x[1]);
	md_free(x0[0]);
	md_free(x0[1]);

	mesh_free(conf);

	EnableWindow(GetDlgItem(hMainDlg, IDC_CALC), TRUE);

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


	GetDlgItemText(hDlg, IDC_COMBO_W2, buf, TEXT_BUF_LENGTH);
	if (lstrcmpi(buf, StrInf) == 0) {
		w2 = -1;
	} else if (0 != StringToDouble(buf, &w2)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	} else {
		if (w2 < 1e-20) {
			mom_error(TEXT("ERROR: The value of `w2` must be > 0"));
			return 1;
		}
		w2 *= 1e-3;
	}


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


	GetDlgItemText(hDlg, IDC_COMBO_PORT, buf, TEXT_BUF_LENGTH);
	if (lstrcmpi(buf, StrAuto) == 0) {
		port_ext = -1;
	} else if (0 != StringToDouble(buf, &port_ext)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	} else {
		port_ext *= 1e-3;
		/* assume w1 is already processed */
		if (port_ext < (w1 * 1e-3)) {
			mom_error(TEXT("small `p` will be supported in the next version"));
			printf("%le\n", port_ext);
			return 1;
		}
	}


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

	GetDlgItemText(hDlg, IDC_COMBO_MESHSTEP, buf, TEXT_BUF_LENGTH);
	if (lstrcmpi(buf, StrAuto) == 0) {
		mesh_step = -1;
	} else if (0 != StringToDouble(buf, &mesh_step)) {
		mom_error(TEXT("unrecognized value:"));
		mom_error(buf);
		return 1;
	} else {
		mesh_step *= 1e-3;
	}

	return 0;
}


INT_PTR CALLBACK AboutWndProc(
	HWND hDlg,
	UINT uMsg,
	WPARAM wParam, LPARAM lParam)
{
	if (uMsg == WM_CLOSE ||
		(uMsg == WM_COMMAND && LOWORD(wParam) == IDOK)) {
		EndDialog(hDlg, 0);
		return TRUE;
	}

	return FALSE;
}


void PaintStaticPicture(HWND hWnd)
{
	HDC hdc, cdc;
	PAINTSTRUCT p;
	RECT dest;
	double wh;
	int w, h;

	hdc = BeginPaint(hWnd, &p);
	cdc = CreateCompatibleDC(hdc);
	SelectObject(cdc, hLayoutBMP);

	GetClientRect(hWnd, &dest);

	wh = (double)LayoutBMP.bmWidth / (double)LayoutBMP.bmHeight;
	w = (int)(wh*(double)dest.bottom);
	if (w < dest.right) {
		h = (int)(((double)w) / wh);
	} else {
		w = dest.right;
		h = (int)(((double)w) / wh);
	}
	w = (int)(0.9*w);
	h = (int)(0.9*h);

	StretchBlt(
		hdc, (dest.right - w)/2, (dest.bottom - h)/2
		, w, h, cdc
		, 0, 0, LayoutBMP.bmWidth, LayoutBMP.bmHeight
		, SRCAND);

	DeleteDC(cdc);
	EndPaint(hWnd, &p);
}

LRESULT CALLBACK StaticImageWndProc(
	HWND hWnd,
	UINT uMsg,
	WPARAM wParam,
	LPARAM lParam)
{
	LRESULT rev = 0;

	switch (uMsg) {
		case WM_PAINT:
			if (GetUpdateRect(hWnd, NULL, FALSE) == 0) {
				break;
			}
			PaintStaticPicture(hWnd);
			break;
		default:
			rev = DefWindowProc(hWnd, uMsg, wParam, lParam);
	}
	return rev;
}


INT_PTR CALLBACK MainWndProc(
	HWND hDlg,
	UINT uMsg,
	WPARAM wParam, LPARAM lParam)
{
	RECT rect;
	switch (uMsg) {
		case WM_INITDIALOG:
			hMainDlg = hDlg;

			SendDlgItemMessage(hDlg, IDC_CHECK_PLOT,
				BM_SETCHECK, BST_CHECKED, 0);

			SendDlgItemMessage(hDlg, IDC_COMBO_W2,
				CB_ADDSTRING, 0, (LPARAM)StrInf);
			SendDlgItemMessage(hDlg, IDC_COMBO_PORT,
				CB_ADDSTRING, 0, (LPARAM)StrAuto);
			SendDlgItemMessage(hDlg, IDC_COMBO_MESHSTEP,
				CB_ADDSTRING, 0, (LPARAM)StrAuto);

			SetDlgItemText(hDlg, IDC_EDIT_W1, TEXT("2.4"));
			SetDlgItemText(hDlg, IDC_COMBO_W2, TEXT("10"));
			SetDlgItemText(hDlg, IDC_EDIT_H, TEXT("0.79"));
			SetDlgItemText(hDlg, IDC_COMBO_PORT, TEXT("2.0"));
			SetDlgItemText(hDlg, IDC_EDIT_EPS, TEXT("2.2"));
			SetDlgItemText(hDlg, IDC_EDIT_D, TEXT("0.0"));
			SetDlgItemText(hDlg, IDC_COMBO_MESHSTEP, StrAuto);

			SetClassLong(hDlg, GCL_HICON,
			    (LONG)LoadIcon(hInst, MAKEINTRESOURCE(IDI_ICON)));
			SetWindowLong(GetDlgItem(hDlg, IDC_IMAGE),
				GWL_WNDPROC, (LONG)StaticImageWndProc);

			GetClientRect(GetDlgItem(hDlg, IDC_LIST_RESULT), &rect);
			SendDlgItemMessage(hDlg, IDC_LIST_RESULT,
				LB_SETHORIZONTALEXTENT, rect.right, 0);

			break;
		case WM_COMMAND:
			switch (LOWORD(wParam)) {
				case IDC_CALC:
					if (Read(hDlg) == 0) {
						hThreadCalc = CreateThread(
							NULL,
							0,
							Calc,
							hDlg,
							0,
							&TIDCalc);
					}
					break;
				case IDC_ABOUT:
					DialogBoxParam(
						hInst,
						MAKEINTRESOURCE(IDD_ABOUT),
						hDlg,
						AboutWndProc,
						0
					);
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

int WINAPI WinMain(HINSTANCE hThis,
		   HINSTANCE hPrevInstance,
		   LPSTR lpCmdLine,
		   int nShowCmd)
{
#ifdef MOM_ENABLE_COMCTL
#ifndef OLD_STYLE_COMCTL
	INITCOMMONCONTROLSEX cc;
#endif
#endif
	hInst = hThis;

	mom_trace = inter_trace;
	mom_error = mom_trace;

	hLayoutBMP = LoadImage(
		hInst, MAKEINTRESOURCE(IDB_LAYOUT),
		IMAGE_BITMAP, 0,0, LR_SHARED | LR_DEFAULTSIZE);

#ifdef MOM_ENABLE_COMCTL
#ifndef OLD_STYLE_COMCTL
	cc.dwSize = sizeof(cc);
	cc.dwICC = ICC_STANDARD_CLASSES;
	InitCommonControlsEx(&cc);
#else
	InitCommonControls();
#endif
#endif

	if (hLayoutBMP == NULL) {
		MessageBox(NULL,
			TEXT("Could not load resource"),
			TEXT("Error"),
			MB_ICONERROR | MB_OK);
		return 1;
	}

	GetObject(hLayoutBMP, sizeof(LayoutBMP), &LayoutBMP);

	return DialogBoxParam(
			hInst,
			MAKEINTRESOURCE(IDD_DLG_MAIN),
			HWND_DESKTOP,
			MainWndProc,
			0
		);
}
