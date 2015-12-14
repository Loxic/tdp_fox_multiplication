#ifndef _STACK_H_
#define _STACK_H_

#include <stdlib.h>

typedef struct stack_node_d{
  void *data;

  struct stack_node_d *next;
} stack_node_t;

stack_node_t * stack_node__create(void *data);
void stack_node__destroy(stack_node_t *s_n);

typedef struct stack_d{
  stack_node_t *top;
  int size;
} stack_t;

stack_t * stack__create();
void stack__destroy(stack_t *s);
void stack__push(void *data, stack_t *s);
void* stack__pop(stack_t *s);
void* stack__peep(stack_t *s);
int stack__is_empty(stack_t *s);


#endif
