# BioRad nsi build script for Windows x86_64 platform
# 
#--------------------------------
!verbose push
!verbose 3

# Maximum compression.
SetCompressor /SOLID lzma


# installation only for current user
!define MULTIUSER_EXECUTIONLEVEL Standard
!define MULTIUSER_INSTALLMODE_INSTDIR "BioRad"
!include MultiUser.nsh

# Define directories.
!define GIT_DIR "."
!define SRC_DIR "${GIT_DIR}\src"
!define BUILD_DIR "${GIT_DIR}\build"
!define DATA_DIR "${GIT_DIR}\database"
!define DOC_DIR "${GIT_DIR}\doc"

!define PYTHON_MAJOR   "3"
!define PYTHON_MINOR   "8"

# The following are derived from the above.
!define PYTHON_VERS    "${PYTHON_MAJOR}.${PYTHON_MINOR}"
!define PYTHON_HK      "Software\Python\PythonCore\${PYTHON_VERS}\InstallPath"
!define PYTHON_HK_64   "Software\Wow6432Node\Python\PythonCore\${PYTHON_VERS}\InstallPath"


# Include the tools we use.
!include MUI2.nsh
!include LogicLib.nsh


Name "BioRad 1.0.5"
Caption "BioRad 1.0.5 Setup"
#InstallDir "$PROGRAMFILES\GeoMop"
OutFile "${GIT_DIR}\dist\BioRad_1.0.5_x86_64.exe"

# Registry key to check for directory (so if you install again, it will 
# overwrite the old one automatically)
InstallDirRegKey HKCU "Software\BioRad" "Install_Dir"

# Request application privileges for Windows Vista and newer
#RequestExecutionLevel admin

#--------------------------------

# Define the different pages.
!define MUI_FINISHPAGE_NOAUTOCLOSE
!insertmacro MUI_PAGE_WELCOME
!insertmacro MUI_PAGE_LICENSE "${GIT_DIR}\LICENSE"
!insertmacro MUI_PAGE_DIRECTORY
!insertmacro MUI_PAGE_INSTFILES
!insertmacro MUI_PAGE_FINISH

!insertmacro MUI_UNPAGE_CONFIRM
!insertmacro MUI_UNPAGE_INSTFILES

# Other settings.
!insertmacro MUI_LANGUAGE "English"

#--------------------------------
# Init

Var PYTHON_EXE
Var PYTHON_SCRIPTS

Function .onInit

  !insertmacro MULTIUSER_INIT

  !define APP_HOME_DIR "$APPDATA\BioRad"

  CheckPython:
    # Check if Python is installed.
    ReadRegStr $PYTHON_EXE HKCU "${PYTHON_HK}" ""

    ${If} $PYTHON_EXE == ""
        ReadRegStr $PYTHON_EXE HKLM "${PYTHON_HK}" ""
    ${Endif}

    # Install Python.
    ${If} $PYTHON_EXE == ""
      MessageBox MB_YESNO|MB_ICONQUESTION "Python ${PYTHON_VERS} 64b is not installed. Do you wish to install it?" IDYES InstallPython
                  Abort
      InstallPython:
        SetOutPath $INSTDIR\prerequisites
        File "${BUILD_DIR}\python-3.8.3-amd64.exe"
        ExecWait 'python-3.8.3-amd64.exe'

        # Check installation.
        Goto CheckPython
    ${Endif}

    # Set the path to the python.exe instead of directory.
    StrCpy $PYTHON_EXE "$PYTHON_EXEpython.exe"

FunctionEnd

Function un.onInit
  !insertmacro MULTIUSER_UNINIT
FunctionEnd

#--------------------------------
# The stuff to install

# These are the programs that are needed by BioRad.
Section "Runtime Environment" SecRuntime
  
  # Section is mandatory.
  SectionIn RO

  # Clean BioRad source, env directories.
  RMDir /r "$INSTDIR\env"

  # Install virtualenv.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\virtualenv-20.0.25-py2.py3-none-any.whl"
  ExecWait '"$PYTHON_EXE" -m pip install "$INSTDIR\prerequisites\virtualenv-20.0.25-py2.py3-none-any.whl"'
  ExecWait '"$PYTHON_EXE" -m virtualenv "$INSTDIR\env"'
  
  SetOutPath $INSTDIR  

  # Copy the src folder.
  File /r "${SRC_DIR}"
  File /r "${DATA_DIR}"
  File /r "${DOC_DIR}"

  # Copy LICENSE, CHANGELOG, VERSION.
  File "${GIT_DIR}\VERSION"
  File "${GIT_DIR}\LICENSE"
  
  # Create folders
  CreateDirectory "$INSTDIR\inputs"
  CreateDirectory "$INSTDIR\outputs"

  # Copy documentation
#  SetOutPath $INSTDIR\doc
#  File "${BUILD_DIR}\Geomop 1.1.0 reference guide.pdf"

  # Set the varible with path to python virtual environment scripts.
  StrCpy $PYTHON_SCRIPTS "$INSTDIR\env\Scripts"

  # Install Matplotlib.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\matplotlib-3.2.2-cp38-cp38-win_amd64.whl"
  ExecWait '"$PYTHON_SCRIPTS\python.exe" -m pip install "$INSTDIR\prerequisites\matplotlib-3.2.2-cp38-cp38-win_amd64.whl"'

  # Install PyYAML.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\PyYAML-5.3.1-cp38-cp38-win_amd64.whl"
  ExecWait '"$PYTHON_SCRIPTS\python.exe" -m pip install "$INSTDIR\prerequisites\PyYAML-5.3.1-cp38-cp38-win_amd64.whl"'

  # Install NumPy.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\numpy-1.19.0-cp38-cp38-win_amd64.whl"
  ExecWait '"$PYTHON_SCRIPTS\python.exe" -m pip install "$INSTDIR\prerequisites\numpy-1.19.0-cp38-cp38-win_amd64.whl"'

  # Install SciPy.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\scipy-1.5.0-cp38-cp38-win_amd64.whl"
  ExecWait '"$PYTHON_SCRIPTS\python.exe" -m pip install "$INSTDIR\prerequisites\scipy-1.5.0-cp38-cp38-win_amd64.whl"'

  # Install pysqlite3.
  SetOutPath $INSTDIR\prerequisites
  File "${BUILD_DIR}\pysqlite3-0.4.3.tar.gz"
  ExecWait '"$PYTHON_SCRIPTS\python.exe" -m pip install "$INSTDIR\prerequisites\pysqlite3-0.4.3.tar.gz"'

SectionEnd




Section "-Batch files" SecBatchFiles

  CreateDirectory "$INSTDIR\bin"
  SetOutPath $INSTDIR\bin

  IfFileExists "$INSTDIR\src\GUI_main.py" 0 +6
    FileOpen $0 "BioRad.bat" w
    FileWrite $0 "@echo off$\r$\n"
    FileWrite $0 'set "PYTHONPATH=$INSTDIR"$\r$\n'
    FileWrite $0 '"$PYTHON_SCRIPTS\python.exe" "$INSTDIR\src\GUI_main.py" %*$\r$\n'
	FileWrite $0 "pause$\r$\n"
    FileClose $0
	
  FileOpen $0 "pythonw.bat" w
  FileWrite $0 "@echo off$\r$\n"
  FileWrite $0 'set "PYTHONPATH=$INSTDIR"$\r$\n'
  FileWrite $0 'start "" "$PYTHON_SCRIPTS\pythonw.exe" %*$\r$\n'
  FileClose $0

SectionEnd


Section "Start Menu shortcuts" SecStartShortcuts

  CreateDirectory "$SMPROGRAMS\BioRad"

  # Make sure this is clean and tidy.
  RMDir /r "$SMPROGRAMS\BioRad"
  CreateDirectory "$SMPROGRAMS\BioRad"
  
  # Uninstall shortcut.
  SetOutPath $INSTDIR
  CreateShortcut "$SMPROGRAMS\BioRad\Uninstall.lnk" "$INSTDIR\uninstall.exe" "" "$INSTDIR\uninstall.exe" 0

  IfFileExists "$INSTDIR\src\GUI_main.py" 0 +3
    SetOutPath $INSTDIR\src
    CreateShortcut "$SMPROGRAMS\BioRad\BioRad.lnk" "$INSTDIR\bin\pythonw.bat" '"$INSTDIR\src\GUI_main.py"' "$INSTDIR\src\BioRad_large.ico" 0

SectionEnd


Section "Desktop icons" SecDesktopIcons

  IfFileExists "$INSTDIR\src\GUI_main.py" 0 +3
    SetOutPath $INSTDIR\src
    CreateShortCut "$DESKTOP\BioRad.lnk" "$INSTDIR\bin\pythonw.bat" '"$INSTDIR\src\GUI_main.py"' "$INSTDIR\src\BioRad_large.ico" 0

SectionEnd


Section /o "Wipe settings" SecWipeSettings

  # Clear all configuration from APPDATA
  RMDir /r "${APP_HOME_DIR}"

SectionEnd


Section "-Default resources data" SecDefaultResourcesData

  # Section is mandatory.
  SectionIn RO

  IfFileExists "${APP_HOME_DIR}" +2 0
    CreateDirectory "${APP_HOME_DIR}"
    # fill data home to default resources data
    #SetOutPath "${APP_HOME_DIR}"
    #File /r "${DATA_DIR}/*"

SectionEnd


Section -post
 
  ; Set output path to the installation directory.
  SetOutPath $INSTDIR
  
  WriteUninstaller "uninstall.exe"

  ; Write the installation path into the registry
  WriteRegStr HKCU SOFTWARE\BioRad "Install_Dir" "$INSTDIR"
  
  ; Write the uninstall keys for Windows
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\BioRad" "DisplayName" "BioRad"
  WriteRegStr HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\BioRad" "UninstallString" "$\"$INSTDIR\uninstall.exe$\""
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\BioRad" "NoModify" 1
  WriteRegDWORD HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\BioRad" "NoRepair" 1
  
SectionEnd


# Uninstaller
Section "Uninstall"
  
  # Remove registry keys
  DeleteRegKey HKCU "Software\Microsoft\Windows\CurrentVersion\Uninstall\BioRad"
  DeleteRegKey HKCU SOFTWARE\BioRad

  # Remove start menu shortcuts.
  RMDir /r "$SMPROGRAMS\BioRad"
  
  # Delete desktop icons.
  Delete "$DESKTOP\BioRad.lnk"

  # Remove BioRad installation directory and all files within.
  RMDir /r "$INSTDIR"

SectionEnd