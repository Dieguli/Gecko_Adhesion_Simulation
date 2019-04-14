//---------------------------------------------------------------------------

#ifndef mainH
#define mainH
//---------------------------------------------------------------------------
#include <Classes.hpp>
#include <Controls.hpp>
#include <StdCtrls.hpp>
#include <Forms.hpp>
#include <ExtCtrls.hpp>
#include <Menus.hpp>
#include <Dialogs.hpp>
#include <gl\gl.h>
#include <gl\glu.h>
#include "View.h"
#include "World.h"

//---------------------------------------------------------------------------
class TmainForm : public TForm
{
__published:	// Von der IDE verwaltete Komponenten
        TPanel *Panel;
        TMainMenu *MainMenu;
        TMenuItem *mHelp;
        TMenuItem *mHelpAbout;
        TTimer *Timer;
        TOpenDialog *OpenDialog;
        TMenuItem *mFile;
        TMenuItem *mFileQuit;
        void __fastcall FormCreate(TObject *Sender);
        void __fastcall FormDestroy(TObject *Sender);
        void __fastcall PanelResize(TObject *Sender);
        void __fastcall FormPaint(TObject *Sender);
        void __fastcall mFileQuitClick(TObject *Sender);
        void __fastcall mHelpAboutClick(TObject *Sender);
        void __fastcall PanelMouseDown(TObject *Sender,
          TMouseButton Button, TShiftState Shift, int X, int Y);
        void __fastcall PanelMouseMove(TObject *Sender, TShiftState Shift,
          int X, int Y);
        void __fastcall PanelMouseUp(TObject *Sender, TMouseButton Button,
          TShiftState Shift, int X, int Y);
        void __fastcall TimerTimer(TObject *Sender);
private:	// Anwender-Deklarationen
  HDC    _hdc;
  HGLRC  _hrc;
  int    _PixelFormat;
  bool   _firstIdle;
  View   _view;
  World *_world;
  Vector3d _worldCenter;
  float    _worldScale;
  bool     _paused;
  bool     _shiftMode;
  int      _frames;
  int      _ticks;
  char     _applicationPath[256];

  void __fastcall OpenGLInit();
  void __fastcall OpenGLDone();
  void __fastcall RenderGLScene();
  void __fastcall IdleLoop(TObject*, bool& done);
  void __fastcall AppMessage(tagMSG &Msg, bool &Handled);

public:		// Anwender-Deklarationen
        __fastcall TmainForm(TComponent* Owner);
};
//---------------------------------------------------------------------------
extern PACKAGE TmainForm *mainForm;
//---------------------------------------------------------------------------
#endif
