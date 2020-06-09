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

set clipboard=unnamedplus

"Comment color
hi Comment ctermfg=Gray

nnoremap <CR> o<Esc>
nnoremap <F2> :set invnumber <CR>
nnoremap <F3> :set invhlsearch <CR>
nnoremap <F4> :set invignorecase <CR>
nnoremap <buffer> <F5> :exec '!python3' shellescape(@%, 1)<CR>

" Commenting blocks of code.
autocmd FileType c,cpp,java,scala let b:comment_leader = '// '
autocmd FileType sh,ruby,python   let b:comment_leader = '# '
autocmd FileType conf,fstab       let b:comment_leader = '# '
autocmd FileType tex              let b:comment_leader = '% '
autocmd FileType mail             let b:comment_leader = '> '
autocmd FileType vim              let b:comment_leader = '" '
autocmd FileType charmm,fortran   let b:comment_leader = '! '
noremap <F9>  :<C-B>silent <C-E>s/^/<C-R>=escape(b:comment_leader,'\/')<CR>/<CR>:nohlsearch<CR>
noremap <F10> :<C-B>silent <C-E>s/^\V<C-R>=escape(b:comment_leader,'\/')<CR>//e<CR>:nohlsearch<CR>
