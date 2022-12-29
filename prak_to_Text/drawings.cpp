#pragma once
#include <SFML/Graphics.hpp>
#include"drawings.h"



int draw_matrix(Eigen::SparseMatrix<int> mat,std::string fn)
{

    int s1 = mat.rows();
    int s2 = mat.cols();
    sf::RenderWindow window(sf::VideoMode(s1, s2), "matrixview");
   
    //sf::RenderWindow window1(sf::VideoMode(s1, s2), "after");
    window.setFramerateLimit(60);
    //window1.setFramerateLimit(60);
    sf::Image image;
    sf::Texture texture;
    sf::Sprite sprite;
    image.create(s1, s2, sf::Color::Black);
    for (int k = 0; k < mat.outerSize(); ++k)
        for (Eigen::SparseMatrix<int>::InnerIterator it(mat, k); it; ++it)
        {
            image.setPixel(it.row(), it.col(), sf::Color::White);
        }
    texture.loadFromImage(image);
    sprite.setTexture(texture);
   


   

    while (window.isOpen() )
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }

        window.clear();
        window.draw(sprite);
        window.display();

    }
    if (fn != "none")
    {
        image.saveToFile(fn + ".jpg");
    }

    return 0;
}
int draw_edges1(std::vector<edge> ed,std::vector<vertex> coord)
{
    int s1 = 720;
    int s2 = 720;
    float midx = s1 / 2;
    float midy = s2 / 2;
    sf::Color color = sf::Color::Blue;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::RenderWindow window(sf::VideoMode(s1, s2), "edges view");
    sf::RenderTexture rendtext;
    rendtext.create(s1, s2);
    window.setFramerateLimit(60);
    //sf::Points.setOrigin(midx, midy);
    std::vector<sf::Vertex> lines;
    lines.resize(2*ed.size());
    for (int k = 0; k< ed.size() - 2; k += 2)
    {
        int c = 10;
        float x1 = float(c * coord[ed[k].start].x + midx);
        float x2 = float(c * coord[ed[k].end].x + midx);
        float y1 = float(c * (-coord[ed[k].start].y) + midy);
        float y2 = float(c * (-coord[ed[k].end].y) + midy);
        sf::Vertex points[] =
        {
            sf::Vertex(sf::Vector2f(x1, y1)),
            sf::Vertex(sf::Vector2f(x2, y2))
            
        };
        points[0].color = color;
        points[1].color = color;
        rendtext.draw(points, 2, sf::Lines);
    }
   // window.clear();
    int k = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear();
        sprite.setTexture(rendtext.getTexture());
        window.draw(sprite);
        window.display();
       // window.clear();
        //window.draw(rendtext);
        //std::vector<sf::Vertex> line;
        //line.resize(2);
        //int k = 0;
        
       // for (int i = 0; i < 2*ed.size()-2; i+=2)
        //{
       
        std::cout << k << std::endl;
        k++;
        //if (k < ed.size())
        //{
        //    window.draw(sprite);
        //    //for(int i=0;i<ed.size)
        //    int c = 10;
        //    float x1 = float(c*coord[ed[k].start].x + midx);
        //    float x2 = float(c*coord[ed[k].end].x + midx);
        //    float y1 = float(c*( - coord[ed[k].start].y) + midy);
        //    float y2 = float(c*( - coord[ed[k].end].y) + midy);
        //    /*float x1 = k;
        //        float x2 = midx;
        //        float y1 = k;
        //        float y2 = midy;*/
        //        //sf::Vertex vertices[2];
        //    sf::Vertex points[] =
        //    {
        //        sf::Vertex(sf::Vector2f(x1, y1)),
        //        sf::Vertex(sf::Vector2f(x2, y2))
        //    };

        //    // vertices[0].position= sf::Vector2f(float(coord[ed[k].start].x - midx), float(coord[ed[k].start].y - midy));
        //     //vertices[1].position= sf::Vector2f(float(coord[ed[k].end].x - midx), float(coord[ed[k].end].y - midy));
        //     //vertices[0].color = sf::Color::Red;
        //     //vertices[1].color = sf::Color::Red;

        //    points[0].color = color;
        //    points[1].color = color;
        //    /*line[0] = sf::Vector2f(coord[ed[k].start].x - midx, coord[ed[k].start].y - midy);
        //    line[1] = sf::Vector2f(coord[ed[k].end].x - midx, coord[ed[k].end].y - midy);*/
        //    window.draw(points, 2, sf::Lines);
        //    window.display();
        //    k++;
        //    texture.create(window.getSize().x, window.getSize().y);
        //    texture.update(window);
        //    sprite.setTexture(texture);
        //   
        //}
       
        //if (k > ed.size())
        //{
        //    window.clear();
        //    window.draw(sprite);
        //    window.display();
        //}
        //}
        

    }
    return 0;

}
int draw_edges(std::vector<edge> ed, std::vector<vertex> coord, std::string fn)
{
    int s1 = 720;
    int s2 = 720;
    float midx = s1 / 2;
    float midy = s2 / 2;
    sf::Color color = sf::Color::Blue;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::CircleShape circ(1.f);
    circ.setFillColor(sf::Color::Red);
    int mass = 5;
    sf::RenderWindow window(sf::VideoMode(s1, s2), "edges view");
    sf::RenderTexture rendtext, rendtext1;
    rendtext.create(s1, s2);
    rendtext1.create(s1, s2);
    window.setFramerateLimit(60);
    //sf::Points.setOrigin(midx, midy);
    std::vector<sf::Vertex> lines;
    lines.resize(2 * ed.size());
    rendtext.clear();
    
    // window.clear();
    int k = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        if (k > ed.size())
        {
            window.close();
        }
        window.clear();
        sprite.setTexture(rendtext.getTexture());
        window.draw(sprite);
        window.display();
        // window.clear();
         //window.draw(rendtext);
         //std::vector<sf::Vertex> line;
         //line.resize(2);
         //int k = 0;

        // for (int i = 0; i < 2*ed.size()-2; i+=2)
         //{

        std::cout << k << std::endl;
       // k++;
        if (k < ed.size())
        {
            window.draw(sprite);
            //for(int i=0;i<ed.size)
            int c = 2000;
           // int c = 10;
            float x1 = float(c * coord[ed[k].start].x + midx-200);
            float x2 = float(c * coord[ed[k].end].x + midx-200);
            float y1 = float(c * (-coord[ed[k].start].y) + midy+200);
            float y2 = float(c * (-coord[ed[k].end].y) + midy+200);

            for (int i = 0; i < mass; i++)
            {
                float xc = i * ((x2 - x1) / mass) + x1;
                float yc = i * ((y2 - y1) / mass) + y1;
                circ.setPosition(xc, yc);
                rendtext.draw(circ);
                //window.draw(circ);
            }
           k++;
            texture.create(window.getSize().x, window.getSize().y);
            texture.update(window);
            sprite.setTexture(texture);
           
        }

        if (k > ed.size())
        {
            window.clear();
            window.draw(sprite);
            window.display();
        }
        


    
    
    }
    if (k < ed.size()) {
        for (int j = 0; j < ed.size(); j++)
        {

        int c = 2000;
        float x1 = float(c * coord[ed[j].start].x + midx - 200);
        float x2 = float(c * coord[ed[j].end].x + midx - 200);
        float y1 = float(c * (-coord[ed[j].start].y) + midy + 200);
        float y2 = float(c * (-coord[ed[j].end].y) + midy + 200);

            for (int i = 0; i < mass; i++)
            {
            float xc = i * ((x2 - x1) / mass) + x1;
            float yc = i * ((y2 - y1) / mass) + y1;
            circ.setPosition(xc, yc);
            rendtext1.draw(circ);
            }


        }
    }
    sprite.setTexture(rendtext1.getTexture());
 
    image=rendtext1.getTexture().copyToImage();
    if (fn != "none")
    {
        image.saveToFile(fn + ".jpg");

    }

    return 0;

}
int draw_coord(std::vector<vertex> coord)
{
    int s1 = 720;
    int s2 = 720;
    float midx = s1 / 2;
    float midy = s2 / 2;
    sf::Color color = sf::Color::Blue;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::CircleShape circ(1.f);
    circ.setFillColor(sf::Color::Red);
    int mass = 5;
    sf::RenderWindow window(sf::VideoMode(s1, s2), "edges view");
    sf::RenderTexture rendtext;
    rendtext.create(s1, s2);
    window.setFramerateLimit(10);
    //sf::Points.setOrigin(midx, midy);
   
    rendtext.clear();
    /*for (int k = 0; k < ed.size() - 2; k += 2)
    {

        int c = 10;
        float x1 = float(c * coord[ed[k].start].x + midx);
        float x2 = float(c * coord[ed[k].end].x + midx);
        float y1 = float(c * (-coord[ed[k].start].y) + midy);
        float y2 = float(c * (-coord[ed[k].end].y) + midy);

        for (int i = 0; i < mass; i++)
        {
            float xc = i*((x2 - x1) / mass) + x1;
            float yc= i*((y2 - y1) / mass) + y1;
            circ.setPosition(xc, yc);
            rendtext.draw(circ);
        }


    }*/
    // window.clear();
    int k = 0;
    while (window.isOpen())
    {
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear();
        sprite.setTexture(rendtext.getTexture());
        window.draw(sprite);
        window.display();
        // window.clear();
         //window.draw(rendtext);
         //std::vector<sf::Vertex> line;
         //line.resize(2);
         //int k = 0;

        // for (int i = 0; i < 2*ed.size()-2; i+=2)
         //{

        std::cout << k << std::endl;
        k++;
        if (k < coord.size())
        {
            window.draw(sprite);
            //for(int i=0;i<ed.size)
            for (int i = 0; i < coord[k].connections.size();i++)
            {
                int c = 10;
                    float x1 = float(c * coord[k].x + midx);
                    float x2 = float(c * coord[coord[k].connections[i]].x + midx);
                    float y1 = float(c * (-coord[k].y) + midy);
                    float y2 = float(c * (-coord[coord[k].connections[i]].y) + midy);

                    for (int j = 0; j < mass; j++)
                    {
                        float xc = j * ((x2 - x1) / mass) + x1;
                            float yc = j * ((y2 - y1) / mass) + y1;
                        circ.setPosition(xc, yc);
                        window.draw(circ);
                    }
                k++;
                texture.create(window.getSize().x, window.getSize().y);
                texture.update(window);
                sprite.setTexture(texture);
            }
        }

        if (k >=coord.size())
        {
            window.clear();
            window.draw(sprite);
            window.display();
        }



    }
    return 0;

}
int draw_graph(std::vector<float> x, std::vector<float> y)
{
    int s1 = 720;
    int s2 = 720;
    float midx = s1 / 2;
    float midy = s2 / 2;
    sf::Color color = sf::Color::Yellow;
    sf::Texture texture;
    sf::Sprite sprite;
    sf::Image image;
    sf::CircleShape circ(1.f);
    circ.setFillColor(color);
    int mass = 5;
    sf::RenderWindow window(sf::VideoMode(s1, s2), "edges view");
    sf::RenderTexture rendtext, rendtext1;
    rendtext.create(s1, s2);
    rendtext1.create(s1, s2);
    window.setFramerateLimit(60);


    while (window.isOpen())
    {
        window.clear();
        sf::Event event;
        while (window.pollEvent(event))
        {
            if (event.type == sf::Event::Closed)
                window.close();
        }
        window.clear();
        sf::Vertex points[] =
             {
                  sf::Vertex(sf::Vector2f(0, midy)),
                  sf::Vertex(sf::Vector2f(s1, midy))
                };
        window.draw(points, 2, sf::Lines);


        for (int i = 0; i < x.size(); i++)
        {
            int c = 3000;
            double cy = 5*pow(10, -7);
            float x1 = float(c * x[i] );
            //float x2 = float(c * coord[ed[j].end].x + midx - 200);
            float y1 = float(cy * (-y[i])+midy); //midy);
            std::cout << y1 <<std:: endl;
            circ.setPosition(x1, y1);
            window.draw(circ);
            //float y2 = float(c * (-coord[ed[j].end].y) + midy + 200);
        }
        
        
        window.display();

    }
    return 0;
}