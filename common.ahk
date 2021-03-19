;--------------------------------------------------------------------
; [Ctrl] + [Alt] shortcuts
;--------------------------------------------------------------------

; [Ctrl] + [Alt] + [D] to open Dictionary
^!D::
    Run, "C:\Users\wangq\Documents\Shortcuts\dict"
    Return

; [Ctrl] + [Alt] + [E] to open Evernote
^!E::
    Run, "C:\Program Files (x86)\Yinxiang Biji\Yinxiang Biji\Evernote.exe"
    Return

; [Ctrl] + [Alt] + [G] to open Google chrome
^!G::
    Run, "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"
    Return

; [Ctrl] + [Alt] + [J] to open JupyterLab in Chrome (makesure JupyterLab is Online)
^!J::
    Run, "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe" "http://localhost:8888/lab"
    Return

; [Ctrl] + [Alt] + [I] to open Inkscape
^!I::
   Run, "C:\Program Files\Inkscape\inkscape.exe"
   Return

; [Ctrl] + [Alt] + [N] to open Marxico to add a New Note
^!N::
    Run, "C:\Program Files (x86)\Maxiang\maxiang.exe"
    Return

; [Ctrl] + [Alt] + [T] to open Terminal
^!T::
    Run, "C:\Users\wangq\Documents\Shortcuts\terminal"
    Return

; [Ctrl] + [Alt] + [V] to open VMD
^!V::
   Run, "C:\Program Files (x86)\University of Illinois\VMD\vmd.exe"
   Return

; [Ctrl] + [Alt] + [W] to open Data Warrior
^!W::
    Run, "C:\Program Files\DataWarrior\DataWarrior.exe"
    Return

;--------------------------------------------------------------------
; Other shortcuts
;--------------------------------------------------------------------

; Use [Pause] key to sleep
Pause::
    DllCall("PowrProf\SetSuspendState", "int", 0, "int", 0, "int", 0)
    Return

; [Win] + [Del] to open trash (recycle bin)
#Del::RUN ::{645ff040-5081-101b-9f08-00aa002f954e} 
; [Win] + [Del] to empty trash (recycle bin)
; #Del::FileRecycleEmpty 

; [Alt] + [W] to close active window
!W::
    WinGetActiveTitle, Title
    WinClose, %Title%
    return

; [Win] + [N] to create a folder named by current date
#N:: 
    Send,+^n %A_YYYY%-%A_MM%-%A_DD%{Enter} 
return
