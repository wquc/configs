"Indentation and other editing options
syntax on
set hlsearch
set tabstop=4
set shiftwidth=4
set expandtab
set number
set showcmd
set showmatch
set incsearch
set foldenable
set t_ut=""
set number
set backspace=indent,eol,start
highlight LineNr term=bold cterm=NONE ctermfg=DarkGrey ctermbg=NONE gui=NONE guifg=DarkGrey guibg=NONE

set cursorline
highlight clear CursorLine
hi CursorLineNr term=bold cterm=bold ctermfg=003 gui=bold

set clipboard=unnamedplus

set background=dark
colorscheme hybrid

"Comment color
"hi Comment ctermfg=Gray

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

call plug#begin()
    Plug 'preservim/nerdtree'
    Plug 'mechatroner/rainbow_csv'
    Plug 'easymotion/vim-easymotion'
    Plug 'vim-airline/vim-airline'
    Plug 'vim-airline/vim-airline-themes'
    Plug 'powerline/powerline'
    Plug 'Yggdroot/indentLine'
call plug#end()

nnoremap <F1> :NERDTreeToggle <CR>
nnoremap <F2> :set invnumber <CR>
nnoremap <F3> :set invhlsearch <CR>
nnoremap <F4> :set invignorecase <CR>

nmap f <Plug>(easymotion-s2)

:let g:airline_section_b = '%{strftime("%H:%M")}'
