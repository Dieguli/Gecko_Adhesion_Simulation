object mainForm: TmainForm
  Left = 0
  Top = 8
  Width = 933
  Height = 726
  Caption = 'Real-Time Physically-Based Environment'
  Color = clBtnFace
  Font.Charset = DEFAULT_CHARSET
  Font.Color = clWindowText
  Font.Height = -11
  Font.Name = 'MS Sans Serif'
  Font.Style = []
  Menu = MainMenu
  OldCreateOrder = False
  Position = poScreenCenter
  Scaled = False
  OnCreate = FormCreate
  OnDestroy = FormDestroy
  OnPaint = FormPaint
  PixelsPerInch = 96
  TextHeight = 13
  object Panel: TPanel
    Left = 0
    Top = 0
    Width = 925
    Height = 672
    Align = alClient
    BevelOuter = bvNone
    FullRepaint = False
    TabOrder = 0
    OnMouseDown = PanelMouseDown
    OnMouseMove = PanelMouseMove
    OnMouseUp = PanelMouseUp
    OnResize = PanelResize
  end
  object MainMenu: TMainMenu
    Left = 17
    Top = 16
    object mFile: TMenuItem
      Caption = '&File'
      object mFileQuit: TMenuItem
        Caption = '&Quit'
        OnClick = mFileQuitClick
      end
    end
    object mHelp: TMenuItem
      Caption = '&Help'
      object mHelpAbout: TMenuItem
        Caption = '&About...'
        ShortCut = 112
        OnClick = mHelpAboutClick
      end
    end
  end
  object Timer: TTimer
    Interval = 10
    OnTimer = TimerTimer
    Left = 73
    Top = 16
  end
  object OpenDialog: TOpenDialog
    Left = 129
    Top = 16
  end
end
