; Use [Pause] key to sleep
Pause::
    DllCall("PowrProf\SetSuspendState", "int", 0, "int", 0, "int", 0)
    Return

; [Ctrl] + [Alt] + [G] to open internet browser
^!G::
    Run, msedge.exe
    Return

; [Ctrl] + [Alt] + [T] to open terminal (need WSL installed)
^!T::
    Run, wsl.exe
    Return

; [Ctrl] + [Alt] + [M] to open Marxico/Maxiang
^!M::
    Run, "C:\Program Files (x86)\Maxiang\maxiang.exe"
    Return

; [Win] + [Del] to open trash (recycle bin)
#Del::RUN ::{645ff040-5081-101b-9f08-00aa002f954e} 

; or [Win] + [Del] to empty trash (recycle bin)
; #Del::FileRecycleEmpty
