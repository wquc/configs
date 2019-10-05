au BufRead,BufNewFile *.inp set filetype=charmm
au! Syntax charmm source $HOME/.vim/syntax/charmm.vim
au BufNewfile,Bufread *.src set filetype=fortran
au BufNewfile,Bufread *.fcm set filetype=fortran
au BufNewfile,Bufread *.str set filetype=charmm
au BufNewfile,Bufread *.namd set filetype=namd

set hlsearch
set tabstop=4
set shiftwidth=4
set expandtab
set number
set showcmd
set showmatch
set incsearch
set foldenable
