# configs
Configure files for commonly used apps.

Notes:
- The `_xxx` file/directory corresponds to the `.xxx` file/directory in `$HOME` folder.
- Always backup the original files !

Although the .rc files are quite self-explanatory, the following will highlight some features of the setup: 
#### BASH configuration
- `ls` : ignore unwanted folders and files, highlight results with color scheme.
- `lll`: list files with human readable file size and sort by creating time.
- `e`  : exit current terminal.
- `p`  : show the current path.
- `..` : goto parent folder.
- `...`: goto grandparent folder.
- `tf` : show the end of an updating file.
- `get`: get the absolute path of a file.
- `CP` : folder copying.
- `RM` : folder removing.
- `cdd`: change directory to a new folder by absolute path.
- `bak`  : make backup for a file
- `bakk` : make backup for a file with date (add %H%M for hours and minutes)
- `fetch`: fetch a PDB file from RCSB site based on its PDB code
- `extract` : extract commonly used compressed files

#### VIM configuration
Besides the syntax highlighting for commonly used programming language. The following hotkeys are enabled:
- `Enter` : Add a new empty line without going to the insert mode and then escape.
- `F2`    : Toggle line number.
- `F3`    : Toggle highlight of search result.
- `F4`    : Toggle case sensitivity of search.
- `F5`    : Run current editing file with python interpreter (default is python3). 
- `F9`    : Comment out current line or selected lines.
- `F10`   : Uncomment current line or selected lines.

#### VMD configuration
After initialization, the script will try to distribute the `main`, `graphics`, `tk console` and `OpenGL display` on the screen with proper position. Besides, the following hotkeys are remapped or enabled:
- `Left Arrow`  : rotate around Y axis by -90.
- `Right Arrow` : rotate around Y axis by  90. 
- `Down Arrow`  : rotate around X axis by  90.
- `Up Arrow`    : rotate around X axis by -90. 
- `Home`        : zoom in by 110%.
- `End`         : zoom out by 90%.
- `Insert`      : Reset current view. Note this is mapped as `=` by default.
- `=`           : Play next frame of the trajectory. Thus pressing `SHIFT` key is no longer required since by default this function is mapped with `+`.
- `Q`           : color top molecule by segment name.
- `W`           : show top molecule as New Cartoon representation.
- `E`           : use AOChalky as the material for top molecule.
- `A`           : turn off axes.
- `D`           : change background color as white.
- `F`           : a combination of all `QWEAD`-mapped operations.
- `F9`          : render current frame using snapshot renderer.
- `F10`         : render current frame using povray renderer. (need povray 3.6 installed)
- `F11`         : render current frame using Tachyon rendered for ambient occlusion effect.
- `F12`         : render current frame using Tachyon rendered for ambient occlusion effect (GPU-accelerated, need properly installed NVIDIA driver).