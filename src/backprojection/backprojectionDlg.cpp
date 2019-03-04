// backprojectionDlg.cpp : implementation file
//

#include "stdafx.h"
#include "backprojection.h"
#include "backprojectionDlg.h"
#include "afxdialogex.h"

#ifdef _DEBUG
#define new DEBUG_NEW
#endif

// CBackProjectionDlg dialog

CBackProjectionDlg::CBackProjectionDlg(CWnd* pParent /*=NULL*/)
    : CSimulationDialog(CBackProjectionDlg::IDD, pParent)
{
    m_hIcon = AfxGetApp()->LoadIcon(IDR_MAINFRAME);
    m_data.params = model::make_default_parameters();
    m_data.stage = model::stage_no;
}

void CBackProjectionDlg::DoDataExchange(CDataExchange* pDX)
{
    CSimulationDialog::DoDataExchange(pDX);
    DDX_Control(pDX, IDC_SRC, m_srcCtrl);
    DDX_Control(pDX, IDC_SRC2, m_dstCtrl);
    DDX_Control(pDX, IDC_SRC3, m_prjCtrl);
    DDX_Text(pDX, IDC_EDIT1, m_data.params.n);
    DDX_Text(pDX, IDC_EDIT2, m_data.params.snr);
}

BEGIN_MESSAGE_MAP(CBackProjectionDlg, CSimulationDialog)
    ON_WM_PAINT()
    ON_WM_QUERYDRAGICON()
    ON_BN_CLICKED(IDC_BUTTON1, &CBackProjectionDlg::OnBnClickedButton1)
    ON_BN_CLICKED(IDC_BUTTON2, &CBackProjectionDlg::OnBnClickedButton2)
END_MESSAGE_MAP()

// CBackProjectionDlg message handlers

BOOL CBackProjectionDlg::OnInitDialog()
{
    CSimulationDialog::OnInitDialog();

    // Set the icon for this dialog.  The framework does this automatically
    //  when the application's main window is not a dialog
    SetIcon(m_hIcon, TRUE);            // Set big icon
    SetIcon(m_hIcon, FALSE);        // Set small icon

    //m_srcCtrl.plot_layer.with(model::make_bmp_plot(m_data.csource));
    m_srcCtrl.plot_layer.with(model::make_bmp_plot(m_data.cnoised));
    m_prjCtrl.plot_layer.with(model::make_bmp_plot(m_data.cprojection));
    m_dstCtrl.plot_layer.with(model::make_bmp_plot(m_data.crecovered));

    // TODO: Add extra initialization here

    return TRUE;  // return TRUE  unless you set the focus to a control
}

// If you add a minimize button to your dialog, you will need the code below
//  to draw the icon.  For MFC applications using the document/view model,
//  this is automatically done for you by the framework.

void CBackProjectionDlg::OnPaint()
{
    if (IsIconic())
    {
        CPaintDC dc(this); // device context for painting

        SendMessage(WM_ICONERASEBKGND, reinterpret_cast<WPARAM>(dc.GetSafeHdc()), 0);

        // Center icon in client rectangle
        int cxIcon = GetSystemMetrics(SM_CXICON);
        int cyIcon = GetSystemMetrics(SM_CYICON);
        CRect rect;
        GetClientRect(&rect);
        int x = (rect.Width() - cxIcon + 1) / 2;
        int y = (rect.Height() - cyIcon + 1) / 2;

        // Draw the icon
        dc.DrawIcon(x, y, m_hIcon);
    }
    else
    {
        CSimulationDialog::OnPaint();
    }
}

// The system calls this function to obtain the cursor to display while the user drags
//  the minimized window.
HCURSOR CBackProjectionDlg::OnQueryDragIcon()
{
    return static_cast<HCURSOR>(m_hIcon);
}


void CBackProjectionDlg::OnBnClickedButton1()
{
    CFileDialog fd(TRUE, TEXT("bmp"));
    if (fd.DoModal() == IDOK)
    {
        CImage img; img.Load(fd.GetPathName());
        auto w = img.GetWidth(), h = img.GetHeight();
        m_data.csource.DeleteObject();
        m_data.csource.Attach(img.Detach());
        m_data.csource.SetBitmapDimension(w, h);
        m_data.source.from_cbitmap(m_data.csource);
        m_data.cnoised.DeleteObject();
        m_srcCtrl.RedrawWindow();

        OnBnClickedButton2();
    }
}


void CBackProjectionDlg::OnBnClickedButton2()
{
    UpdateData(TRUE);

    m_data.noised = m_data.source;

    double m = 0, e = 0, en = 0;
    for (size_t i = 0; i < m_data.source.h; ++i)
    for (size_t j = 0; j < m_data.source.w; ++j)
        m += m_data.source.data[i][j];
    m /= m_data.source.h * m_data.source.w;
    for (size_t i = 0; i < m_data.source.h; ++i)
    for (size_t j = 0; j < m_data.source.w; ++j)
    {
        double n = 0;
        for (size_t _ = 0; _ < 12; ++_) n += rand() / (RAND_MAX + 1.) - 0.5;
        m_data.noised.data[i][j] = n / 12;
        double d = m_data.source.data[i][j] - m;
        e += d * d;
        en += n * n;
    }

    double a = std::sqrt(m_data.params.snr / 100. * e / en);

    for (size_t i = 0; i < m_data.source.h; ++i)
    for (size_t j = 0; j < m_data.source.w; ++j)
    {
        m_data.noised.data[i][j] =
            std::abs(m_data.source.data[i][j] + a * m_data.noised.data[i][j]);
    }

    m_data.noised.to_cbitmap(m_data.cnoised, 1, false);
    m_srcCtrl.RedrawWindow();

    model::backprojection bp(m_data.params);
    bp.project(m_data.noised, m_data.projection);
    m_data.projection.to_cbitmap(m_data.cprojection, 1, false);
    m_prjCtrl.RedrawWindow();
    bp.backproject(m_data.projection, m_data.recovered);
    m_data.recovered.to_cbitmap(m_data.crecovered, 1, false);
    m_dstCtrl.RedrawWindow();
}
