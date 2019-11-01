"Syntax highlighting
au BufRead,BufNewFile *.inp set filetype=charmm
au! Syntax charmm source $HOME/.vim/syntax/charmm.vim
au BufNewfile,Bufread *.src set filetype=fortran
au BufNewfile,Bufread *.fcm set filetype=fortran
au BufNewfile,Bufread *.str set filetype=charmm
au BufNewfile,Bufread *.namd set filetype=namd

"Indentation and other editing options
set hlsearch
set tabstop=4
set shiftwidth=4
set expandtab
set number
set showcmd
set showmatch
set incsearch
set foldenable

set number
highlight LineNr term=bold cterm=NONE ctermfg=DarkGrey ctermbg=NONE gui=NONE guifg=DarkGrey guibg=NONE

set cursorline
highlight CursorLine cterm=NONE ctermbg=NONE ctermfg=NONE guibg=NONE guifg=NONE

"se mouse+=a
"se mouse=a
set clipboard=unnamedplus

map <F2> :set invnumber <CR>
map <F3> :set invhlsearch <CR>
map <F4> :set invignorecase <CR>
