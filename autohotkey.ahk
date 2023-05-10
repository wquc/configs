;--------------------------------------------------------------------
; Open Applications
;--------------------------------------------------------------------
RunOrActivate(ApplicationPath, ApplicationAhkExe)
{
    Process, Exist, %ApplicationAhkExe%
    If Not ErrorLevel ; errorlevel will = 0 if process doesn't exist
        Run, %ApplicationPath%
    Else
        WinActivate, ahk_exe %ApplicationAhkExe%
    Return
}

; [Ctrl] + [Alt] + [G] to open Google chrome
^!G::
    RunOrActivate("C:\Program Files\Google\Chrome\Application\chrome.exe", "chrome.exe")
    Return

; [Ctrl] + [Alt] + [D] to open Dictionary
^!D::
    Run, "C:\Users\wq\Documents\Shortcuts\Dictionary"
    Return

;--------------------------------------------------------------------
; Utilities
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
