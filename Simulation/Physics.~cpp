//---------------------------------------------------------------------------

#include <vcl.h>
#pragma hdrstop
USERES("Physics.res");
USEFORM("main.cpp", mainForm);
USEUNIT("View.cpp");
USEUNIT("Vertex.cpp");
USEUNIT("Vector3d.cpp");
USEUNIT("Matrix3d.cpp");
USEUNIT("Tetra.cpp");
USEUNIT("TetraMesh.cpp");
USEUNIT("World.cpp");
USEUNIT("Draw.cpp");
USEUNIT("Bounds3d.cpp");
USEUNIT("TetraArray.cpp");
USEUNIT("PhysicsParams.cpp");
USEUNIT("VertexArray.cpp");
USEUNIT("Utils.cpp");
USEFORM("DlgAbout.cpp", AboutDlg);
//---------------------------------------------------------------------------
WINAPI WinMain(HINSTANCE, HINSTANCE, LPSTR, int)
{
        try
        {
                 Application->Initialize();
                 Application->CreateForm(__classid(TmainForm), &mainForm);
                 Application->CreateForm(__classid(TAboutDlg), &AboutDlg);
                 Application->Run();
        }
        catch (Exception &exception)
        {
                 Application->ShowException(&exception);
        }
        return 0;
}
//---------------------------------------------------------------------------
