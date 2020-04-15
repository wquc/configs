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
