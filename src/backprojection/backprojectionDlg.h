// backprojectionDlg.h : header file
//

#include <util/common/gui/SimulationDialog.h>
#include <util/common/gui/PlotControl.h>

#include "model.h"

#pragma once

// CBackProjectionDlg dialog
class CBackProjectionDlg : public CSimulationDialog
{
// Construction
public:
    CBackProjectionDlg(CWnd* pParent = NULL);    // standard constructor

// Dialog Data
    enum { IDD = IDD_BACKPROJECTION_DIALOG };

    protected:
    virtual void DoDataExchange(CDataExchange* pDX);    // DDX/DDV support

// Implementation
protected:
    HICON m_hIcon;

    // Generated message map functions
    virtual BOOL OnInitDialog();
    afx_msg void OnPaint();
    afx_msg HCURSOR OnQueryDragIcon();
    DECLARE_MESSAGE_MAP()
public:
    CPlotControl m_srcCtrl;
    CPlotControl m_prjCtrl;
    CPlotControl m_dstCtrl;
    model::model_data m_data;
    afx_msg void OnBnClickedButton1();
    afx_msg void OnBnClickedButton2();
    enum { algo_backproj = 0, algo_kaczmarz, algo_maxentropy };
    int m_nAlgo;
};
