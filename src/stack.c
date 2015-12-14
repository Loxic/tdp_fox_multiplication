#include "stack.h"

stack_node_t * stack_node__create(void *data){
  stack_node_t *r = malloc(sizeof(*r));
  r->data = data;
  r->next = NULL;
  
  return r;
}

void stack_node__destroy( stack_node_t *s_n){
  free(s_n);
}

stack_t * stack__create(){
  stack_t *r = malloc(sizeof(*r));
  r->top = NULL;
  r->size = 0;

  return r;
}

void stack__destroy(stack_t *s){
  stack_node_t *n = s->top;
  while(n){
    stack_node_t *tmp = n;
    n = n->next;
    stack_node__destroy(tmp);
  }
  free(s);
}

void stack__push(void *data,  stack_t *s){
  stack_node_t *n = stack_node__create(data);
  n->next=s->top;
  s->top=n;
  s->size++;
}

void * stack__pop(stack_t *s){
  void *r = s->top->data;
  stack_node_t *tmp = s->top; 
  s->top=s->top->next;
  s->size--;
  free(tmp);
  return r;
}

void * stack__peep(stack_t *s){
  return s->top->data;
}

int stack__is_empty(stack_t *s){
  return !(s->size);
}
