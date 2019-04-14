//---------------------------------------------------------------------------

#include <vcl.h>
#include <sysutils.hpp>
#include "main.h"
#include "draw.h"
#include "physicsParams.h"
#include "DlgAbout.h"
#include "Bounds3d.h"

#define NEAR_CLIP_PLANE 50
#define FAR_CLIP_PLANE 5000

//---------------------------------------------------------------------------
#pragma package(smart_init)
#pragma resource "*.dfm"
TmainForm *mainForm;

//---------------------------------------------------------------------------
void __fastcall TmainForm::OpenGLInit()//---------------------------------------------------------------------------
{
  _hdc = GetDC(Panel->Handle);
  PIXELFORMATDESCRIPTOR pfd = {
    sizeof(PIXELFORMATDESCRIPTOR),
    1,
    PFD_DRAW_TO_WINDOW | PFD_SUPPORT_OPENGL | PFD_DOUBLEBUFFER,
    PFD_TYPE_RGBA,
    24,
    0,0,0,0,0,0,
    0,0,
    0,0,0,0,0,
    32,
    0,
    0,
    PFD_MAIN_PLANE,
    0,
    0,0,
  };
  _PixelFormat = ChoosePixelFormat(_hdc, &pfd);  SetPixelFormat(_hdc, _PixelFormat, &pfd);
  _hrc = wglCreateContext(_hdc);
  if(_hrc == NULL)
    ShowMessage("OpenGL problem: Error initializing Form object");
  if(wglMakeCurrent(_hdc, _hrc) == false)
    ShowMessage("OpenGL problem: Error initializing Form object");

  Set8087CW(0x133f);  // Disable floating point exceptions                      // recommended when using OpenGL}

//---------------------------------------------------------------------------
void __fastcall TmainForm::OpenGLDone()//---------------------------------------------------------------------------
{
  wglMakeCurrent(_hdc, NULL);
  wglDeleteContext(_hrc);
}


//---------------------------------------------------------------------------
__fastcall TmainForm::TmainForm(TComponent* Owner)
        : TForm(Owner)
//---------------------------------------------------------------------------
{
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::FormCreate(TObject *Sender)
//---------------------------------------------------------------------------
{
  OpenGLInit();
  draw3d.initOpenGL();
  
  PanelResize(NULL);
  physicsParamsInit();

  strcpy(_applicationPath, ExtractFilePath(Application->ExeName).c_str());

  AnsiString s; s.printf("%s\\Models", _applicationPath);
  OpenDialog->InitialDir = s;

  _world = new World(_applicationPath);
  Bounds3d bounds; _world->getBounds(bounds);

  _world->getBounds(bounds);
  if (bounds.isEmpty()) {
    _worldCenter.set(0,0.2,0);
    _worldScale = 0.2 * Height;
  }
  else {
    _worldCenter.add(bounds._min, bounds._max); _worldCenter.scale(0.5);
    _worldScale = Height * 0.2 / bounds._min.dist(bounds._max);
  }
  _view.Reset();

  _paused = true;
  _shiftMode = false;
  _frames = 0;
  _ticks = GetTickCount();

  _firstIdle = true;
  Application->OnIdle = IdleLoop;
  Application->OnMessage = AppMessage;
}


//---------------------------------------------------------------------------
void __fastcall TmainForm::FormDestroy(TObject *Sender)
//---------------------------------------------------------------------------
{
  OpenGLDone();
  delete _world;
}


//---------------------------------------------------------------------------
void __fastcall TmainForm::PanelResize(TObject *Sender)
//---------------------------------------------------------------------------
{
  // left, bottom, width, height
  glViewport(0, 0, Panel->Width, Panel->Height);

  glMatrixMode(GL_PROJECTION);
  glLoadIdentity();

  gluPerspective(40.0, (GLfloat)(Panel->Width)/(GLfloat)(Panel->Height), NEAR_CLIP_PLANE, FAR_CLIP_PLANE);
}


//---------------------------------------------------------------------------
void __fastcall TmainForm::RenderGLScene()
//---------------------------------------------------------------------------
{
  glClear(GL_COLOR_BUFFER_BIT | GL_DEPTH_BUFFER_BIT);

  _view.LookAt(_worldScale, _worldCenter);

  _world->draw();
  glFlush();
  SwapBuffers(_hdc);
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::FormPaint(TObject *Sender)
//---------------------------------------------------------------------------
{
  RenderGLScene();
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::IdleLoop(TObject*, bool& done)
//---------------------------------------------------------------------------
{
  if (_firstIdle) {
    _firstIdle = false;
    RenderGLScene();
  }
}


//---------------------------------------------------------------------------
void __fastcall TmainForm::AppMessage(tagMSG &Msg, bool &Handled)
//---------------------------------------------------------------------------
{
  if (Msg.message == WM_KEYDOWN)
  {
    if (Msg.wParam == 'P') {
      _paused = !_paused;
      Handled = true;
    }
  }
}



//----- Main Menu -----------------------------------------------------


//---------------------------------------------------------------------------
void __fastcall TmainForm::mFileQuitClick(TObject *Sender)
//---------------------------------------------------------------------------
{
  Application->Terminate();
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::mHelpAboutClick(TObject *Sender)
//---------------------------------------------------------------------------
{
  AboutDlg->ShowModal();
}

// ----------------- Mouse handling ------------------------------------------

//---------------------------------------------------------------------------
void __fastcall TmainForm::PanelMouseDown(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
//---------------------------------------------------------------------------
{
  _shiftMode = Shift.Contains(ssShift);

  if (_shiftMode) {
    if (_paused) _paused = false;
    _world->mouseDown(X,Y);
  }
  else
    _view.Click(X,Y);
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::PanelMouseMove(TObject *Sender,
      TShiftState Shift, int X, int Y)
//---------------------------------------------------------------------------
{
  if (_shiftMode) {
    _world->mouseMove(X,Y);
  }
  else {
    bool pressed = true;
    if (Shift.Contains(ssLeft))
      _view.DragRotate(X,Y);
    else if (Shift.Contains(ssMiddle))
      _view.DragTranslate(X,Y);
    else if (Shift.Contains(ssRight))
      _view.DragScale(X,Y);
    else pressed = false;
    if (pressed) FormPaint(NULL);
  }
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::PanelMouseUp(TObject *Sender,
      TMouseButton Button, TShiftState Shift, int X, int Y)
//---------------------------------------------------------------------------
{
  if (_shiftMode)
    _world->mouseUp(X,Y);
}

//---------------------------------------------------------------------------
void __fastcall TmainForm::TimerTimer(TObject *Sender)
//---------------------------------------------------------------------------
{
  if (!_paused)
    _world->animate();
  FormPaint(NULL);
}

