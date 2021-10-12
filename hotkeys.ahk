;--------------------------------------------------------------------
; [Ctrl] + [Alt] shortcuts
;--------------------------------------------------------------------

; [Ctrl] + [Alt] + [G] to open Google chrome
^!G::
    Run, "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe"
    Return

; [Ctrl] + [Alt] + [J] to open JupyterLab in Chrome (makesure JupyterLab is Online)
^!J::
    Run, "C:\Program Files (x86)\Google\Chrome\Application\chrome.exe" "http://localhost:8888/lab"
    Return

; [Ctrl] + [Alt] + [N] to open Marxico to add a New Note
^!N::
    Run, "C:\Program Files (x86)\Maxiang\maxiang.exe"
    Return

; [Ctrl] + [Alt] + [T] to open Terminal
^!T::
    Run, "C:\Users\wangq\Documents\Shortcuts\terminal"
    Return

;--------------------------------------------------------------------
; Other shortcuts
;--------------------------------------------------------------------

; [Alt] + [Q] to close active window
!Q::
    WinGetActiveTitle, Title
    WinClose, %Title%
    return

; [Win] + [N] to create a folder named by current date
#N:: 
    Send,+^n %A_YYYY%-%A_MM%-%A_DD%{Enter} 
return
