//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop

#include "DlgAbout.h"
//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TAboutDlg *AboutDlg;
//---------------------------------------------------------------------------
__fastcall TAboutDlg::TAboutDlg(TComponent* Owner)
        : TForm(Owner)
{
}
//---------------------------------------------------------------------------
void __fastcall TAboutDlg::LabelClick(TObject *Sender)
{
  TLabel *label = (TLabel *)Sender;
  ShellExecute(0,NULL,label->Caption.c_str(),NULL,NULL,SW_SHOW);
}
//---------------------------------------------------------------------------

