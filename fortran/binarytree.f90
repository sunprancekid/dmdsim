! AUTHOR : Matthew A Dorsey
! DATE : 2022-12-20
! FILENAME : binarytree.f90
! PURPOSE : contains binary tree efficiency technique module and methods

module binarytree
implicit none

type :: node 
    integer :: rnode, lnode, pnode
end type node 

contains

! method that initializes a node type to zero pointers
function reset_node () result (n)
    implicit none 
    type(node) :: n
    n%rnode = 0 
    n%lnode = 0
    n%pnode = 0
end function reset_node

! method that initializes an array used for binary tree method
! (all pointers in a binary tree are initialized to zero)
subroutine initialize_binarytree(tree, rootnode)
    implicit none
    type(node), dimension(:), intent(out) :: tree ! event tree containing pointers
    integer, intent(out) :: rootnode ! initial node of binary tree
    integer :: i ! indexing parameter

    ! reset initial node of binary tree 
    rootnode = 0
    ! loop through each node and reset pointers to zero 
    do i = 1, size(tree)
        tree(i) = reset_node()
    end do 
end subroutine initialize_binarytree

! finds the soonest event in a binary event tree
! (soonest event is left most node)
! stored with binary tree module
integer function findnextevent(tree, rootnode)
    implicit none
    type(node), dimension(:), intent(inout) :: tree ! binary tree
    integer, intent(in) :: rootnode ! inital pointer used for binary tree
    integer :: nextnode

    ! start the search for the minimum time with the root node 
    nextnode = rootnode
    do
        if (tree(nextnode)%lnode == 0) exit
        nextnode = tree(nextnode)%lnode
    enddo 
    findnextevent = nextnode
end function findnextevent

! subroutine that removes a branch from the event tree
! algorithm ensures that event order doesn't change in event tree
! adapted from Smith et al. MD for Polymeric Fluids
! stored with binary tree module
subroutine delbranch(tree, rootnode, nonode)
    implicit none
    type(node), dimension(:), intent(inout) :: tree ! binary tree
    integer, intent(inout) :: rootnode ! initial pointer in binary tree
    integer, intent(in) :: nonode ! mold id of particle whose node is being deleted from tree
    integer :: ns, np ! pointers used for relinking

    ! determine the relationship of the deleted node to other nodes in the tree structure
    if (tree(nonode)%rnode == 0) then
        ! CASE I: the deleted node is followed on the right branch by a null event
        ! SOL: pnode of nonode should be linked to lnode of nonode (note lnode of nonode can be
        ! either a null event or another branch)
        ns = tree(nonode)%lnode
        np = tree(nonode)%pnode
        if (ns /= 0) tree(ns)%pnode = np
        if (np /= 0) then 
            if (tree(np)%lnode == nonode) then 
                tree(np)%lnode = ns
            else
                tree(np)%rnode = ns 
            endif
        endif
    else if (tree(nonode)%lnode == 0) then ! 
        ! CASE II: the deleted node contains a null event on the left branch 
        ! and a non-null event on the right branch
        ! SOL: pnode of nonode should be linked to rnode of nonode
        ns = tree(nonode)%rnode
        np = tree(nonode)%pnode
        tree(ns)%pnode = np
        if (np /= 0) then 
            if (tree(np)%lnode == nonode) then 
                tree(np)%lnode = ns
            else
                tree(np)%rnode = ns 
            endif
        endif
    else if (tree(tree(nonode)%rnode)%lnode == 0) then
        ! CASE III: the deleted node contains non-null events on the left and right branches
        ! while the right branch contains a null left branch, which indicates that the right branch
        ! of the right branch is the smallest event time
        ! SOL: Since the event time on the right of nonode is larger than the event time on the left,
        ! rnode of nonode is designated as the successor. The null event of lnode of rnode of nonode
        ! is replaced with lnode of nonode      ! link pnode of nonode to rnode to nonode
        ns = tree(nonode)%rnode
        ! link the left branch of nonode to the left branch of the successor node (null event)
        np = tree(nonode)%lnode 
        tree(ns)%lnode = np 
        tree(np)%pnode = ns 
        ! link the successor node to the pointer node of nonode
        np = tree(nonode)%pnode
        tree(ns)%pnode = np
        if (np /= 0) then 
            if (tree(np)%lnode == nonode) then 
                tree(np)%lnode = ns
            else
                tree(np)%rnode = ns 
            endif
        endif
    else 
        ! CASE IV: last case, most generic solution required. The right branch of nonode has
        ! a non-null left branch
        ! SOL: search left branch of right branch of nonode for the minimum event 
        ! find the node whose event is closest to nonode
        ns = tree(tree(nonode)%rnode)%lnode
        do
            if (tree(ns)%lnode == 0) exit
            ns = tree(ns)%lnode
        enddo 
        ! replace the successor node with the successor node's right branch
        np = tree(ns)%pnode 
        tree(np)%lnode = tree(ns)%rnode 
        if (tree(ns)%rnode /= 0) tree(tree(ns)%rnode)%pnode = np 
        ! replace nonode with the successor node
        ! link the right branch of nonode to the right branch of the successor node
        np = tree(nonode)%rnode 
        tree(np)%pnode = ns 
        tree(ns)%rnode = np
        ! link the left branch of nonode to the left branch of the successor node (null event)
        np = tree(nonode)%lnode 
        tree(np)%pnode = ns 
        tree(ns)%lnode = np 
        ! link the successor node to nonode's predossesor node 
        np = tree(nonode)%pnode
        tree(ns)%pnode = np
        if (np /= 0) then 
            if (tree(np)%lnode == nonode) then 
                tree(np)%lnode = ns
            else
                tree(np)%rnode = ns 
            endif
        endif
    endif

    if (nonode == rootnode) then ! if nonode is the rootnode
        ! reset the rootnode as the successor node
        rootnode = ns 
    endif

    ! reset nonode
    tree(nonode) = reset_node()
end subroutine delbranch
    
end module binarytree